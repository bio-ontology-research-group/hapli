"""Tests that schema works in pydantic-less environments.

We discovered this when running on a cluster env (unimatrix01's esmfold-env)
that had torch + fair-esm but not pydantic. The original fallback assigned
`BaseModel = object`, which doesn't accept kwargs in __init__, so calling
`EpistasisResidual(transcript=..., haplotype=..., ...)` crashed at runtime.

This test exercises the fallback path by running the schema module in a
subprocess with pydantic blocked from import, ensuring it still produces
usable model instances and a working `model_dump()`.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap


_SCRIPT = textwrap.dedent("""
    import builtins, sys
    orig_import = builtins.__import__
    def block(name, *a, **k):
        if name == 'pydantic' or name.startswith('pydantic.'):
            raise ImportError('blocked for test')
        return orig_import(name, *a, **k)
    builtins.__import__ = block

    # Force re-import in case any module pre-imported pydantic.
    for m in list(sys.modules):
        if m == 'pydantic' or m.startswith('pydantic.') or m.startswith('hapli.core.schema'):
            del sys.modules[m]

    from hapli.core import schema

    assert schema._HAS_PYDANTIC is False, "expected pydantic-less mode"

    # Schema models should accept kwargs and store them as attributes.
    r = schema.EpistasisResidual(
        transcript='T1', haplotype=1, n_variants=2,
        s_additive=-58.5, s_joint=-1.4, residual=57.1,
        flagged=True,
    )
    assert r.transcript == 'T1'
    assert r.haplotype == 1
    assert r.n_variants == 2
    assert r.s_additive == -58.5
    assert r.s_joint == -1.4
    assert r.residual == 57.1
    assert r.flagged is True

    # model_dump() must produce a JSON-serialisable dict.
    d = r.model_dump()
    assert d['transcript'] == 'T1'
    assert d['s_additive'] == -58.5

    # Nested models also work.
    p = schema.PresenceCall(status='intact', source='liftoff', coverage=0.95)
    assert p.status == 'intact'
    assert p.coverage == 0.95
    pd = p.model_dump()
    assert pd['status'] == 'intact'

    # AnalysisResult round-trip.
    ev = schema.GeneEvidence(gene='G1')
    ev.presence = {'hap1': p}
    ev.epistasis = [r]
    a = schema.AnalysisResult(gene='G1', sample='S1', region='chr1:1-100',
                              transcripts=[], evidence=ev)
    payload = a.model_dump()
    assert payload['gene'] == 'G1'
    # The nested PresenceCall and EpistasisResidual must have been recursively dumped.
    pres = payload['evidence']['presence']['hap1']
    assert pres['status'] == 'intact', f"nested dump broke: {pres}"
    epi = payload['evidence']['epistasis'][0]
    assert epi['residual'] == 57.1

    print('OK')
""")


def test_schema_works_without_pydantic():
    """Run the schema module in a subprocess with pydantic blocked."""
    res = subprocess.run(
        [sys.executable, "-c", _SCRIPT],
        capture_output=True, text=True, timeout=30,
    )
    assert res.returncode == 0, (
        f"pydantic-less schema fallback broke:\nstdout:\n{res.stdout}\nstderr:\n{res.stderr}"
    )
    assert "OK" in res.stdout
