import json
import os
import urllib.request
import logging
from typing import Dict, Any

class LLMInterpreter:
    def __init__(self, api_key: str = None, model: str = "openai/gpt-3.5-turbo"):
        self.api_key = api_key or os.environ.get("OPENROUTER_API_KEY")
        self.model = model
        self.logger = logging.getLogger(__name__)

    def interpret(self, analysis_data: Dict[str, Any]) -> str:
        if not self.api_key:
            self.logger.warning("No OpenRouter API key provided. Skipping LLM interpretation.")
            return "Interpretation skipped (No API Key)."

        prompt = self._construct_prompt(analysis_data)
        
        try:
            response = self._call_openrouter(prompt)
            return response
        except Exception as e:
            self.logger.error(f"LLM Interpretation failed: {e}")
            return f"Error during interpretation: {e}"

    def _construct_prompt(self, data: Dict[str, Any]) -> str:
        gene = data.get('gene', 'Unknown')
        sample = data.get('sample', 'Unknown')
        transcripts = data.get('transcripts', [])
        ev = data.get('evidence') or {}

        lines = [f"Analyze the functional impact of variants in gene {gene} for sample {sample}."]

        # ── Per-haplotype presence (Liftoff) ───────────────────────────────
        presence = ev.get('presence') or {}
        if presence:
            lines.append("\nGene presence per haplotype (from Liftoff):")
            for hap in ("hap1", "hap2"):
                pc = presence.get(hap)
                if isinstance(pc, dict):
                    status = pc.get('status', '?')
                    sid = pc.get('sequence_identity')
                    extra = f", lift identity {sid:.2f}" if isinstance(sid, (int, float)) else ""
                    lines.append(f"  - {hap}: {status}{extra}")

        # ── Diploid summary ────────────────────────────────────────────────
        dip = ev.get('diploid') or {}
        if dip:
            lines.append("\nDiploid summary:")
            for k in ("hap1_score", "hap2_score", "min_score", "max_score"):
                v = dip.get(k)
                if v is not None:
                    lines.append(f"  - {k}: {v:.3f}" if isinstance(v, (int, float)) else f"  - {k}: {v}")
            if dip.get('compound_het_lof'):
                lines.append("  - compound_het_lof: TRUE (both haplotypes independently loss-of-function)")

        # ── Per-variant consequences (bcftools/csq) ────────────────────────
        cons = ev.get('consequence') or []
        if cons:
            lines.append(f"\nPer-variant consequences ({len(cons)} records):")
            for c in cons:
                hap = c.get('haplotype', '?')
                kind = c.get('consequence', '?')
                aa = c.get('amino_acid_change') or ''
                clnsig = c.get('clnsig') or ''
                parts = [f"hap{hap}", kind]
                if aa: parts.append(aa)
                if clnsig: parts.append(f"ClinVar: {clnsig}")
                lines.append("  - " + ", ".join(parts))

        # ── Per-haplotype protein diff ─────────────────────────────────────
        prot = ev.get('protein') or []
        if prot:
            lines.append("\nProtein-level diff per haplotype:")
            for p in prot:
                hap = p.get('haplotype', '?')
                ident = p.get('identity')
                pts = p.get('premature_stop_at')
                fsr = p.get('frameshift_region') or {}
                fs = "frameshift restored" if fsr.get('restored') else "no frameshift"
                parts = [f"hap{hap}"]
                if isinstance(ident, (int, float)): parts.append(f"identity {ident:.3f}")
                if pts is not None: parts.append(f"premature stop at AA {pts}")
                parts.append(fs)
                lines.append("  - " + ", ".join(parts))

        # ── Epistasis residual (hapli's novelty) ───────────────────────────
        epi = ev.get('epistasis') or []
        if epi:
            lines.append("\nEpistasis residual (additive − joint protein-language-model scores):")
            for r in epi:
                hap = r.get('haplotype', '?')
                resid = r.get('residual')
                flagged = r.get('flagged')
                if isinstance(resid, (int, float)):
                    tag = " [FLAGGED]" if flagged else ""
                    lines.append(
                        f"  - hap{hap}: residual {resid:+.2f} "
                        f"(S_add={r.get('s_additive'):+.2f}, S_joint={r.get('s_joint'):+.2f}){tag}"
                    )

        # ── Constraints (gnomAD / ClinGen) ─────────────────────────────────
        constr = ev.get('constraints') or {}
        if constr and any(constr.values()):
            lines.append("\nGene constraint priors:")
            for k in ("pli", "mis_z", "oe_lof", "clingen_haploinsufficiency"):
                v = constr.get(k)
                if v is not None:
                    lines.append(f"  - {k}: {v}")

        # ── Fallback to schema-v1 alignment block if no evidence present ──
        if not ev and transcripts:
            lines.append("\nAlignment results (legacy schema):")
            for t in transcripts:
                tid = t.get('id', '?')
                lines.append(f"Transcript {tid}:")
                for hap, info in (t.get('alignments') or {}).items():
                    if isinstance(info, str):
                        lines.append(f"  - {hap}: {info}")
                    else:
                        identity = info.get('identity', 0)
                        lines.append(f"  - {hap}: Identity {identity:.2%}")

        lines.append(
            "\nGiven the above, explain the functional impact of this gene on each "
            "haplotype, whether it's loss-of-function on one or both alleles, and "
            "any epistatic signals (e.g., large |residual| implies variants "
            "interact non-additively). Be concise and actionable for a clinician."
        )
        return "\n".join(lines)

    def _call_openrouter(self, prompt: str) -> str:
        url = "https://openrouter.ai/api/v1/chat/completions"
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
            "HTTP-Referer": "https://hapli.ai", # Optional
        }
        
        payload = {
            "model": self.model,
            "messages": [
                {"role": "system", "content": "You are a bioinformatics expert specializing in variant interpretation."},
                {"role": "user", "content": prompt}
            ]
        }
        
        req = urllib.request.Request(
            url, 
            data=json.dumps(payload).encode('utf-8'), 
            headers=headers, 
            method="POST"
        )
        
        with urllib.request.urlopen(req) as response:
            result = json.loads(response.read().decode('utf-8'))
            return result['choices'][0]['message']['content']
