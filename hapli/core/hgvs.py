"""
HGVS protein-level variant parser + haplotype-protein constructor.

Covers the variant classes that show up in MAVE scoresets and in haplotype
construction from VCF + GFF3 pipelines:

  p.Arg200Trp              missense
  p.Arg200=                synonymous
  p.Arg200Ter  /  p.R200*  nonsense (stop-gain)
  p.Arg200fs               frameshift (anchor only, effect unspecified)
  p.Arg200GlyfsTer10       frameshift with alt first residue + premature stop
  p.Arg200del              single-residue in-frame deletion
  p.Arg200_Lys210del       range in-frame deletion
  p.Met1_Val50del          N-terminal range deletion (effectively whole-gene
                           LoF for short proteins)
  p.Arg200delinsGlyLeu     single-residue in-frame delins
  p.[Arg200Trp;Lys210Glu]  multi-variant haplotype (any mix of the above)

Applying a list of variants to a reference protein yields a synthetic
"haplotype protein" the downstream scorer can treat uniformly, whether the
underlying variants are SNV missense, VCF indels, or symbolic SVs that have
been pre-processed into HGVS form by Phase 0–1.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable


class VariantKind(str, Enum):
    MISSENSE = "missense"
    SYNONYMOUS = "synonymous"
    NONSENSE = "nonsense"
    FRAMESHIFT = "frameshift"
    DEL = "del"                  # in-frame deletion (single or range)
    DELINS = "delins"
    UNKNOWN = "unknown"


# Three-letter <-> one-letter amino acids (HGVS is three-letter by default)
_AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O",
    "Ter": "*", "*": "*",
}


def _aa_to_1(code: str) -> str | None:
    """Accept either 3-letter or 1-letter amino acid, return 1-letter, else None."""
    if not code:
        return None
    if code in _AA3_TO_1:
        return _AA3_TO_1[code]
    if len(code) == 1 and code in "ACDEFGHIKLMNPQRSTVWYU*":
        return code
    return None


# ─────────────────────────────────────────────────────────────────────────────
# One variant — structured result of parsing
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class ProteinVariant:
    """One HGVS-p variant. Coordinates are 1-based inclusive.

    For missense/synonymous/nonsense: start == end, ref_aa is 1-letter, alt_aa is 1-letter.
    For range DEL: start..end span, ref_aa/alt_aa unused.
    For frameshift: start is the anchor, `fs_first_alt_aa` is the first new aa (may be None),
                    `fs_ter_offset` is the +N offset to the first in-frame stop (may be None).
    For delins: start..end span plus `alt_seq` (inserted aa string).
    """

    kind: VariantKind
    start: int                              # 1-based inclusive
    end: int                                # 1-based inclusive
    raw: str                                # original HGVS string
    ref_aa: str | None = None               # single aa (missense/nonsense)
    alt_aa: str | None = None               # single aa (missense/nonsense)
    alt_seq: str | None = None              # multi-aa insertion (delins)
    fs_first_alt_aa: str | None = None      # frameshift
    fs_ter_offset: int | None = None        # frameshift (+N to PTC)


# ─────────────────────────────────────────────────────────────────────────────
# Regex catalogue
# ─────────────────────────────────────────────────────────────────────────────
_THREE_OR_ONE = r"(?:[A-Z][a-z]{2}|[A-Z])"
_POS = r"\d+"

_RE_MISSENSE = re.compile(
    rf"^({_THREE_OR_ONE})({_POS})({_THREE_OR_ONE}|\*|Ter)$"
)
_RE_SYNONYMOUS = re.compile(
    rf"^({_THREE_OR_ONE})({_POS})=$"
)
_RE_NONSENSE = re.compile(
    rf"^({_THREE_OR_ONE})({_POS})(?:Ter|\*)$"
)
_RE_DEL_SINGLE = re.compile(
    rf"^({_THREE_OR_ONE})({_POS})del$"
)
_RE_DEL_RANGE = re.compile(
    rf"^({_THREE_OR_ONE})({_POS})_({_THREE_OR_ONE})({_POS})del$"
)
_RE_DELINS_SINGLE = re.compile(
    rf"^({_THREE_OR_ONE})({_POS})delins([A-Za-z]+)$"
)
_RE_DELINS_RANGE = re.compile(
    rf"^({_THREE_OR_ONE})({_POS})_({_THREE_OR_ONE})({_POS})delins([A-Za-z]+)$"
)
_RE_FRAMESHIFT = re.compile(
    rf"^({_THREE_OR_ONE})({_POS})(?:({_THREE_OR_ONE})?fs(?:Ter|\*)?({_POS})?)$"
)


def _parse_single(raw: str) -> ProteinVariant | None:
    """Parse one variant token without the leading `p.` prefix."""
    raw = raw.strip()

    # Order matters: nonsense matches missense regex too (alt='Ter'), so test nonsense first.
    if (m := _RE_NONSENSE.match(raw)):
        ref3, pos = m.groups()
        ref1 = _aa_to_1(ref3)
        if ref1 is None:
            return None
        pos_i = int(pos)
        return ProteinVariant(
            kind=VariantKind.NONSENSE, start=pos_i, end=pos_i,
            raw=raw, ref_aa=ref1, alt_aa="*",
        )

    if (m := _RE_SYNONYMOUS.match(raw)):
        ref3, pos = m.groups()
        ref1 = _aa_to_1(ref3)
        if ref1 is None:
            return None
        pos_i = int(pos)
        return ProteinVariant(
            kind=VariantKind.SYNONYMOUS, start=pos_i, end=pos_i,
            raw=raw, ref_aa=ref1, alt_aa=ref1,
        )

    if (m := _RE_DEL_RANGE.match(raw)):
        ref3a, pa, ref3b, pb = m.groups()
        a, b = int(pa), int(pb)
        if b < a:
            return None
        return ProteinVariant(
            kind=VariantKind.DEL, start=a, end=b, raw=raw,
            ref_aa=_aa_to_1(ref3a), alt_aa=None,
        )

    if (m := _RE_DEL_SINGLE.match(raw)):
        ref3, pos = m.groups()
        pos_i = int(pos)
        return ProteinVariant(
            kind=VariantKind.DEL, start=pos_i, end=pos_i, raw=raw,
            ref_aa=_aa_to_1(ref3), alt_aa=None,
        )

    if (m := _RE_DELINS_RANGE.match(raw)):
        ref3a, pa, ref3b, pb, alt_text = m.groups()
        a, b = int(pa), int(pb)
        if b < a:
            return None
        alt_seq = _decode_aa_string(alt_text)
        if alt_seq is None:
            return None
        return ProteinVariant(
            kind=VariantKind.DELINS, start=a, end=b, raw=raw,
            alt_seq=alt_seq,
        )

    if (m := _RE_DELINS_SINGLE.match(raw)):
        ref3, pos, alt_text = m.groups()
        pos_i = int(pos)
        alt_seq = _decode_aa_string(alt_text)
        if alt_seq is None:
            return None
        return ProteinVariant(
            kind=VariantKind.DELINS, start=pos_i, end=pos_i, raw=raw,
            alt_seq=alt_seq,
        )

    if (m := _RE_FRAMESHIFT.match(raw)):
        ref3, pos, fs_alt3, ter_offset = m.groups()
        pos_i = int(pos)
        return ProteinVariant(
            kind=VariantKind.FRAMESHIFT, start=pos_i, end=pos_i, raw=raw,
            ref_aa=_aa_to_1(ref3),
            fs_first_alt_aa=_aa_to_1(fs_alt3) if fs_alt3 else None,
            fs_ter_offset=int(ter_offset) if ter_offset else None,
        )

    if (m := _RE_MISSENSE.match(raw)):
        ref3, pos, alt3 = m.groups()
        ref1 = _aa_to_1(ref3)
        alt1 = _aa_to_1(alt3)
        if ref1 is None or alt1 is None:
            return None
        pos_i = int(pos)
        if ref1 == alt1:
            return ProteinVariant(
                kind=VariantKind.SYNONYMOUS, start=pos_i, end=pos_i,
                raw=raw, ref_aa=ref1, alt_aa=ref1,
            )
        return ProteinVariant(
            kind=VariantKind.MISSENSE, start=pos_i, end=pos_i,
            raw=raw, ref_aa=ref1, alt_aa=alt1,
        )

    return None


def _decode_aa_string(s: str) -> str | None:
    """Turn a string of three-letter codes or one-letter codes into one-letter."""
    if not s:
        return None
    # All-uppercase 1-letter: simplest case.
    if all(c in "ACDEFGHIKLMNPQRSTVWYU*" for c in s):
        return s
    # Three-letter repeated codes (sizeof 3 each).
    if len(s) % 3 != 0:
        return None
    out = []
    for i in range(0, len(s), 3):
        one = _aa_to_1(s[i : i + 3])
        if one is None:
            return None
        out.append(one)
    return "".join(out)


def parse_hgvs_pro(notation: str) -> list[ProteinVariant]:
    """Parse one HGVS protein notation string. Returns possibly-empty list of variants.

    Accepts bracketed multi-variant form `p.[X;Y;Z]` and the bare `p.X` form.
    Returns an empty list if nothing could be parsed.
    """
    if not notation or notation == "NA":
        return []
    s = notation.strip()
    if s.startswith("p."):
        s = s[2:]
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1]
    if s.startswith("[") and s.endswith("]"):
        tokens = [t.strip() for t in s[1:-1].split(";") if t.strip()]
    else:
        tokens = [s]
    out: list[ProteinVariant] = []
    for t in tokens:
        v = _parse_single(t)
        if v is not None:
            out.append(v)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Apply variants to a reference protein → haplotype protein
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class HaplotypeProtein:
    """The constructed haplotype protein + provenance flags."""

    seq: str                                 # translated haplotype protein (may be empty / truncated)
    variants: list[ProteinVariant] = field(default_factory=list)
    categories: list[str] = field(default_factory=list)
    applied: list[ProteinVariant] = field(default_factory=list)
    skipped: list[ProteinVariant] = field(default_factory=list)
    premature_stop_at: int | None = None      # 1-based aa index within hap_seq
    length_changed: bool = False
    overall_category: str = "unmodified"


def apply_variants(ref_seq: str, variants: Iterable[ProteinVariant]) -> HaplotypeProtein:
    """Apply the variants to the reference protein in reverse-coordinate order
    (so an upstream frameshift/deletion doesn't shift positions of downstream
    variants). Returns the constructed haplotype protein plus diagnostics.

    Any variant whose coordinates don't match the reference is skipped with a
    note, not raised — real MAVE data contains rows referencing annotations
    from older UniProt releases.
    """
    vlist = list(variants)
    ordered = sorted(vlist, key=lambda v: -v.start)      # reverse coord order

    skipped: list[ProteinVariant] = []
    applied: list[ProteinVariant] = []
    categories: list[str] = []
    hap = list(ref_seq)
    ptc: int | None = None

    def reject(v: ProteinVariant, why: str) -> None:
        skipped.append(v)
        categories.append(f"skipped:{v.kind.value}:{why}")

    for v in ordered:
        if v.start < 1 or v.end > len(ref_seq):
            reject(v, "out_of_range"); continue

        if v.kind == VariantKind.SYNONYMOUS:
            # no-op
            applied.append(v); categories.append("synonymous"); continue

        if v.kind == VariantKind.MISSENSE:
            if hap[v.start - 1] != v.ref_aa:
                reject(v, "ref_mismatch"); continue
            hap[v.start - 1] = v.alt_aa or "X"
            applied.append(v); categories.append("missense"); continue

        if v.kind == VariantKind.NONSENSE:
            if hap[v.start - 1] != v.ref_aa:
                reject(v, "ref_mismatch"); continue
            # Truncate the haplotype at the stop position (keep nothing after).
            hap = hap[: v.start - 1]
            ptc = v.start
            applied.append(v); categories.append("nonsense"); continue

        if v.kind == VariantKind.DEL:
            # In-frame deletion of range [start, end] (1-based inclusive).
            if v.ref_aa is not None and hap[v.start - 1] != v.ref_aa:
                reject(v, "ref_mismatch"); continue
            del hap[v.start - 1 : v.end]
            applied.append(v); categories.append("del"); continue

        if v.kind == VariantKind.DELINS:
            alt = v.alt_seq or ""
            # Replace range [start, end] with alt_seq.
            hap[v.start - 1 : v.end] = list(alt)
            applied.append(v); categories.append("delins"); continue

        if v.kind == VariantKind.FRAMESHIFT:
            # Model as: truncate at the anchor, append the altered first aa (if
            # given), then append `fs_ter_offset - 1` unknown residues (or
            # simply truncate if ter_offset is unknown). This is a conservative
            # structural approximation — the real downstream sequence could be
            # anything, but we don't have the DNA. For scoring purposes the
            # key signals are (a) length change (captured), (b) premature stop
            # (captured), (c) divergent sequence in the tail (not captured —
            # left as 'X' stand-ins which ESM will score as low-probability).
            new_tail: list[str] = []
            if v.fs_first_alt_aa:
                new_tail.append(v.fs_first_alt_aa)
            if v.fs_ter_offset is not None and v.fs_ter_offset > 1:
                new_tail.extend(["X"] * (v.fs_ter_offset - 1))
            hap = hap[: v.start - 1] + new_tail
            ptc = v.start + (v.fs_ter_offset or 0)
            applied.append(v); categories.append("frameshift"); continue

        reject(v, "unsupported_kind")

    hap_seq = "".join(hap)
    length_changed = len(hap_seq) != len(ref_seq)

    # Derive an overall category for routing downstream.
    unique = set(categories)
    applied_cats = {c for c in unique if not c.startswith("skipped:")}
    if not applied_cats:
        overall = "unmodified"
    elif "nonsense" in applied_cats or "frameshift" in applied_cats:
        overall = "lof"
    elif "del" in applied_cats or "delins" in applied_cats:
        overall = "indel"
    elif len(applied_cats) == 1 and "missense" in applied_cats:
        overall = "missense_single" if len([v for v in applied if v.kind == VariantKind.MISSENSE]) == 1 else "missense_multi"
    elif applied_cats == {"missense"}:
        overall = "missense_multi"
    elif applied_cats == {"synonymous"}:
        overall = "synonymous"
    else:
        overall = "mixed"

    return HaplotypeProtein(
        seq=hap_seq,
        variants=vlist,
        categories=categories,
        applied=applied,
        skipped=skipped,
        premature_stop_at=ptc,
        length_changed=length_changed,
        overall_category=overall,
    )


__all__ = [
    "HaplotypeProtein",
    "ProteinVariant",
    "VariantKind",
    "apply_variants",
    "parse_hgvs_pro",
]
