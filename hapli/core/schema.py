"""
Schema v2 for the hapli per-gene evidence bundle.

The output JSON retains the schema-v1 top-level keys (`gene`, `sample`,
`region`, `transcripts`) so the existing TUI and LLM consumers keep working,
and nests all new information under a top-level `evidence` key. Each stage
of the pipeline writes into its own sub-bundle:

    evidence.presence         — Liftoff / TOGA per-haplotype gene-presence call
    evidence.consequence      — bcftools csq per-haplotype per-transcript HGVS
    evidence.splice           — SpliceAI Δ scores per variant per junction
    evidence.lof              — LOFTEE HC/LC + NMD-escape classification
    evidence.missense_agg     — AlphaMissense summed over the gene, per hap
    evidence.protein          — reference vs. haplotype protein diff + epistasis residual
    evidence.diploid          — two-number diploid summary + constraint priors

Consumers read defensively via dict.get() (or the helper `read_evidence`)
because any given stage may be absent (opt-in, not run, or not applicable
to this case).

Pydantic v2 is used for writers that want validation; dataclass-style
access works too via `.model_dump()`. The contract is the JSON shape,
not the Python class — consumers that cannot take a pydantic dependency
can treat it as plain dict-of-dicts.
"""

from __future__ import annotations

from typing import Any, Literal, Optional

try:
    from pydantic import BaseModel, Field, ConfigDict
    _HAS_PYDANTIC = True
except ImportError:                                  # pragma: no cover
    _HAS_PYDANTIC = False

    class _ShimBase:
        """Fallback when pydantic isn't installed.

        Accepts any keyword args and exposes class-level defaults as instance
        attributes. Defines `model_dump()` so consumers that pickle the
        evidence bundle still work. The contract is the JSON shape — losing
        runtime validation is acceptable on minimal envs (e.g. cluster scratch
        environments where torch + esm are present but the project's full
        dep set isn't).
        """

        def __init__(self, **kwargs: Any) -> None:
            for k, v in kwargs.items():
                setattr(self, k, v)

        def model_dump(self) -> dict:
            out: dict = {}
            for k in dir(self):
                if k.startswith("_"):
                    continue
                v = getattr(self, k, None)
                if callable(v):
                    continue
                if hasattr(v, "model_dump"):
                    v = v.model_dump()
                elif isinstance(v, list):
                    v = [x.model_dump() if hasattr(x, "model_dump") else x for x in v]
                elif isinstance(v, dict):
                    v = {kk: (vv.model_dump() if hasattr(vv, "model_dump") else vv)
                         for kk, vv in v.items()}
                out[k] = v
            return out

    BaseModel = _ShimBase  # type: ignore
    Field = lambda *a, **k: (k.get("default_factory") and k["default_factory"]()) or k.get("default")  # type: ignore
    ConfigDict = lambda **k: None  # type: ignore


SCHEMA_VERSION = "2.0"


# ─────────────────────────────────────────────────────────────────────────────
# Presence (Liftoff / TOGA) — per-haplotype per-gene
# ─────────────────────────────────────────────────────────────────────────────
PresenceStatus = Literal[
    "intact",            # lifted cleanly, coverage≈1, matches_ref_protein=True
    "partial",           # lifted but coverage <1 or exons missing
    "low_identity",      # lifted but sequence ID below threshold
    "duplicated",        # extra_copy_number > 0
    "deleted",           # unmapped — gene absent from haplotype
    "uncertain",         # lifted but ambiguous (valid_ORF=False etc.)
    "not_run",           # step skipped
]


class PresenceCall(BaseModel):  # type: ignore
    """Per-haplotype Liftoff/TOGA result for one gene."""
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    status: PresenceStatus = "not_run"
    source: str = "liftoff"                          # tool name
    # Coordinates on the *haplotype* (1-based inclusive). None if absent.
    seqid: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    strand: Optional[str] = None
    # Liftoff-provided QC. None if unknown.
    coverage: Optional[float] = None
    sequence_identity: Optional[float] = None
    valid_orfs: Optional[int] = None
    matches_ref_protein: Optional[bool] = None
    extra_copy_number: Optional[int] = None
    # Free-form flags forwarded by the tool (partial_mapping=True, low_identity=True, …)
    flags: dict[str, Any] = Field(default_factory=dict) if _HAS_PYDANTIC else {}
    # Per-transcript presence/ORF-validity, keyed by transcript ID.
    transcripts: dict[str, dict[str, Any]] = Field(default_factory=dict) if _HAS_PYDANTIC else {}


# ─────────────────────────────────────────────────────────────────────────────
# Consequence (bcftools csq) — per-variant per-haplotype per-transcript
# ─────────────────────────────────────────────────────────────────────────────
class ConsequenceCall(BaseModel):  # type: ignore
    """One BCSQ record from bcftools csq."""
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""
    haplotype: int = 0                                # 1 or 2
    transcript: str = ""
    consequence: str = ""                             # e.g. "missense", "stop_gained"
    hgvs_c: Optional[str] = None                      # "c.100A>G"
    hgvs_p: Optional[str] = None                      # "p.Lys34Arg"
    # Joint-variant compound annotation when bcftools csq reports a combined effect
    compound_id: Optional[str] = None


# ─────────────────────────────────────────────────────────────────────────────
# Splice / LoF / Missense-aggregate — lightweight containers
# ─────────────────────────────────────────────────────────────────────────────
class SpliceCall(BaseModel):  # type: ignore
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""
    gene: str = ""
    ds_ag: float = 0.0                                # Δ acceptor gain
    ds_al: float = 0.0                                # Δ acceptor loss
    ds_dg: float = 0.0                                # Δ donor gain
    ds_dl: float = 0.0                                # Δ donor loss
    max_delta: float = 0.0


class LoFCall(BaseModel):  # type: ignore
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""
    haplotype: int = 0
    transcript: str = ""
    confidence: Literal["HC", "LC", "not_lof", "not_run"] = "not_run"
    filter_reasons: list[str] = Field(default_factory=list) if _HAS_PYDANTIC else []
    flags: list[str] = Field(default_factory=list) if _HAS_PYDANTIC else []


class MissenseAggregate(BaseModel):  # type: ignore
    """Per-haplotype per-transcript AlphaMissense aggregate."""
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    transcript: str = ""
    haplotype: int = 0
    n_missense: int = 0
    max_am: float = 0.0
    mean_am: float = 0.0
    sum_am: float = 0.0
    source: str = "alphamissense_v1"


# ─────────────────────────────────────────────────────────────────────────────
# Protein — reference-vs-haplotype diff + epistasis residual
# ─────────────────────────────────────────────────────────────────────────────
class ProteinDiff(BaseModel):  # type: ignore
    """Per-haplotype translated-protein comparison against the reference."""
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    transcript: str = ""
    haplotype: int = 0
    ref_length: int = 0
    hap_length: int = 0
    identity: float = 1.0
    # Structured per-residue diff: list of (ref_pos, ref_aa, hap_pos, hap_aa) tuples.
    substitutions: list[dict[str, Any]] = Field(default_factory=list) if _HAS_PYDANTIC else []
    premature_stop_at: Optional[int] = None           # aa position of PTC, if any
    frameshift_region: Optional[dict[str, Any]] = None  # {start, end, restored}


class EpistasisResidual(BaseModel):  # type: ignore
    """S_joint − S_additive for ≥2 coding variants on the same haplotype."""
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    transcript: str = ""
    haplotype: int = 0
    n_variants: int = 0
    s_additive: float = 0.0
    s_joint: float = 0.0
    residual: float = 0.0
    additive_source: str = "alphamissense"            # or "esm1v_wildtype_marginals"
    joint_source: str = "esm2_pseudo_likelihood"      # or "saprot"
    flagged: bool = False                             # |residual| > threshold


# ─────────────────────────────────────────────────────────────────────────────
# Diploid — two-number summary + constraint priors (no dominance model)
# ─────────────────────────────────────────────────────────────────────────────
class ConstraintPriors(BaseModel):  # type: ignore
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    pli: Optional[float] = None
    mis_z: Optional[float] = None
    oe_lof: Optional[float] = None
    oe_mis: Optional[float] = None
    clingen_haploinsufficient: Optional[bool] = None
    clingen_triplosensitive: Optional[bool] = None


class DiploidReport(BaseModel):  # type: ignore
    """Deliberately NOT a single call — two numbers + context for downstream."""
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    hap1_score: Optional[float] = None
    hap2_score: Optional[float] = None
    min_score: Optional[float] = None
    max_score: Optional[float] = None
    # True iff both haplotypes independently called LoF/absent (compound-het LoF pattern)
    compound_het_lof: bool = False
    constraints: ConstraintPriors = Field(default_factory=ConstraintPriors) if _HAS_PYDANTIC else None  # type: ignore


# ─────────────────────────────────────────────────────────────────────────────
# The full evidence bundle
# ─────────────────────────────────────────────────────────────────────────────
class GeneEvidence(BaseModel):  # type: ignore
    """Everything hapli knows about one gene in one sample."""
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    gene: str = ""
    # Per-haplotype presence call (keys: 'hap1', 'hap2')
    presence: dict[str, PresenceCall] = Field(default_factory=dict) if _HAS_PYDANTIC else {}
    # Per-transcript per-haplotype records
    consequence: list[ConsequenceCall] = Field(default_factory=list) if _HAS_PYDANTIC else []
    splice: list[SpliceCall] = Field(default_factory=list) if _HAS_PYDANTIC else []
    lof: list[LoFCall] = Field(default_factory=list) if _HAS_PYDANTIC else []
    missense_agg: list[MissenseAggregate] = Field(default_factory=list) if _HAS_PYDANTIC else []
    protein: list[ProteinDiff] = Field(default_factory=list) if _HAS_PYDANTIC else []
    epistasis: list[EpistasisResidual] = Field(default_factory=list) if _HAS_PYDANTIC else []
    diploid: DiploidReport = Field(default_factory=DiploidReport) if _HAS_PYDANTIC else None  # type: ignore


# ─────────────────────────────────────────────────────────────────────────────
# Top-level analysis output — schema v1 fields + evidence v2
# ─────────────────────────────────────────────────────────────────────────────
class AnalysisResult(BaseModel):  # type: ignore
    """The JSON emitted by `hapli analyze` / `hapli assess`.

    v1 top-level fields (`gene`, `sample`, `region`, `transcripts`) are kept
    exactly to preserve the TUI + LLM consumers; everything new lives under
    `evidence`.
    """
    if _HAS_PYDANTIC:
        model_config = ConfigDict(extra="allow")

    schema_version: str = SCHEMA_VERSION
    gene: str = ""
    sample: str = ""
    region: str = ""
    transcripts: list[dict[str, Any]] = Field(default_factory=list) if _HAS_PYDANTIC else []
    evidence: GeneEvidence = Field(default_factory=GeneEvidence) if _HAS_PYDANTIC else None  # type: ignore


# ─────────────────────────────────────────────────────────────────────────────
# Helpers for defensive reading
# ─────────────────────────────────────────────────────────────────────────────
def read_evidence(payload: dict[str, Any], *path: str, default: Any = None) -> Any:
    """Walk a dotted path through the dict-form of an AnalysisResult.

    Example:
        read_evidence(json_payload, 'evidence', 'presence', 'hap1', 'status')
    Returns `default` at the first missing key so callers don't need nested ifs.
    """
    cur: Any = payload
    for key in path:
        if cur is None:
            return default
        if isinstance(cur, dict):
            cur = cur.get(key, default)
        else:
            return default
    return cur if cur is not None else default
