"""
ESM2 scoring primitives for the epistasis-residual computation.

Exposes two operations over protein sequences:

  * `per_position_log_probs(seq)` — return (len_seq,) numpy array where entry i
    is log p(seq[i] | seq_{≠i}), computed by feeding the unmasked sequence to
    the forward pass and reading the log-softmax at position i for the token
    that was actually there. This is the "wildtype-marginals" pseudo-
    log-likelihood per residue used by Meier 2021 (ESM-1v) and the subsequent
    ESM2 variant-effect literature.

  * `alt_log_probs(ref_seq, positions, alt_residues)` — for each (pos, alt)
    pair, return log p(alt | ref_seq_{≠pos}). Used to score a variant in a
    *reference* context (the additive baseline in the epistasis residual).

Results are cached on disk under $HAPLI_CACHE_DIR/esm (or ~/.cache/hapli/esm)
keyed by (model name, sha256(sequence)). The cache is fast to hit and tiny
(one .npy per sequence), so a population scan over haplotypes that share
most of their proteome is essentially free after the first pass.

Torch and fair-esm are *optional* dependencies (the `ml` extra). If they are
not importable, `ESM2Scorer` cannot be constructed and callers must handle
`EsmNotAvailable`.
"""

from __future__ import annotations

import hashlib
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


class EsmNotAvailable(ImportError):
    """Raised when the `ml` extra (torch + fair-esm) is not installed."""


def _default_cache_dir() -> Path:
    root = os.environ.get("HAPLI_CACHE_DIR") or str(Path.home() / ".cache" / "hapli")
    return Path(root) / "esm"


def _sha256(s: str) -> str:
    return hashlib.sha256(s.encode("ascii", errors="replace")).hexdigest()


@dataclass
class ESM2Scorer:
    """Per-position log-probability scorer backed by a fair-esm ESM2 model.

    Parameters
    ----------
    model_name : one of the fair-esm pretrained names. Defaults to the tiny
                 `esm2_t6_8M_UR50D` (7.5M params, ~30 MB checkpoint) which is
                 small enough to run unit tests on CPU. For real benchmarks
                 use `esm2_t33_650M_UR50D` or larger with GPU.
    device     : 'cpu' or 'cuda'. Defaults to cuda if available.
    cache_dir  : directory for per-sequence log-probability tensors.
    """

    model_name: str = "esm2_t6_8M_UR50D"
    device: str = ""                                  # auto-detected if empty
    cache_dir: Path = Path()                          # auto-defaulted if empty

    def __post_init__(self) -> None:
        self._logger = logging.getLogger(__name__)
        if not self.cache_dir or self.cache_dir == Path():
            self.cache_dir = _default_cache_dir()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        try:
            import torch
        except ImportError as exc:                    # pragma: no cover
            raise EsmNotAvailable(
                "`torch` is required for ESM2 scoring. Install the `ml` extra: "
                "`uv sync --extra ml` or `pip install hapli[ml]`."
            ) from exc
        try:
            import esm
        except ImportError as exc:                    # pragma: no cover
            raise EsmNotAvailable(
                "`fair-esm` is required. Install the `ml` extra."
            ) from exc

        self._torch = torch
        self._esm = esm

        if not self.device:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"

        loader = getattr(esm.pretrained, self.model_name, None)
        if loader is None:
            raise ValueError(f"unknown ESM model {self.model_name!r}")
        model, alphabet = loader()
        model.eval()
        model = model.to(self.device)
        self._model = model
        self._alphabet = alphabet
        self._batch_converter = alphabet.get_batch_converter()
        self._logger.info(
            "ESM2 model %s loaded on %s (%.1f M params)",
            self.model_name, self.device,
            sum(p.numel() for p in model.parameters()) / 1e6,
        )

    # ─────────────────────────────────────────────────────────────────────
    # Internal helpers
    # ─────────────────────────────────────────────────────────────────────
    def _token_idx(self, aa: str) -> int:
        idx = self._alphabet.get_idx(aa)
        if idx == self._alphabet.unk_idx:
            raise ValueError(f"unknown residue {aa!r}")
        return idx

    def _cache_path(self, key: str, kind: str = "lp") -> Path:
        return self.cache_dir / f"{self.model_name}.{kind}.{key}.npy"

    # ─────────────────────────────────────────────────────────────────────
    # Public API
    # ─────────────────────────────────────────────────────────────────────
    def per_position_log_probs(self, seq: str) -> np.ndarray:
        """Return shape (len(seq),) of log p(seq[i] | seq) via a full-sequence
        forward pass. Results are cached by (model, sha256(seq)).
        """
        if not seq:
            return np.zeros(0, dtype=np.float32)

        cache_key = _sha256(seq)
        cache_path = self._cache_path(cache_key, "lp")
        if cache_path.exists():
            try:
                return np.load(cache_path)
            except Exception:                         # pragma: no cover
                pass

        lp = self._forward_log_probs(seq)[
            np.arange(len(seq)),
            [self._token_idx(aa) for aa in seq],
        ]
        np.save(cache_path, lp.astype(np.float32))
        return lp

    def pseudo_log_likelihood(self, seq: str) -> float:
        return float(self.per_position_log_probs(seq).sum())

    def alt_log_probs(
        self,
        seq: str,
        positions: Sequence[int],                     # 1-based positions in seq
        alts: Sequence[str],
    ) -> np.ndarray:
        """log p(alt_k | seq) for each (position, alt) pair, using the full
        (unmutated) sequence as context. One forward pass shared across all
        positions. Returns (len(positions),) log-probs.
        """
        if not positions:
            return np.zeros(0, dtype=np.float32)
        full = self._forward_log_probs(seq)
        out = np.zeros(len(positions), dtype=np.float32)
        for k, (pos, alt) in enumerate(zip(positions, alts)):
            out[k] = full[pos - 1, self._token_idx(alt)]
        return out

    # ─────────────────────────────────────────────────────────────────────
    # Core forward pass — returns (L, V) log-softmax over the vocabulary for
    # the L residues of `seq` (BOS/EOS trimmed).
    # ─────────────────────────────────────────────────────────────────────
    def _forward_log_probs(self, seq: str) -> np.ndarray:
        torch = self._torch
        _, _, toks = self._batch_converter([("query", seq)])
        toks = toks.to(self.device)
        with torch.no_grad():
            out = self._model(toks, repr_layers=[], return_contacts=False)
        logits = out["logits"][0]                    # (L+2, V) incl. BOS/EOS
        log_softmax = torch.log_softmax(logits, dim=-1)
        # Drop BOS (idx 0) and EOS (idx -1) rows.
        return log_softmax[1:-1].detach().cpu().numpy()


def sequences_differ_only_at(
    ref: str, hap: str, positions: Iterable[int]
) -> bool:
    """Sanity check: return True iff ref and hap agree at every position EXCEPT
    the (1-based) `positions`. Used by the epistasis module to verify that its
    additive/joint decomposition only touches substitution positions.
    """
    if len(ref) != len(hap):
        return False
    positions = set(positions)
    for i in range(1, len(ref) + 1):
        if i in positions:
            continue
        if ref[i - 1] != hap[i - 1]:
            return False
    return True
