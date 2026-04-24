from textual.app import App, ComposeResult
from textual.widgets import Header, Footer, Tree, Static, Label
from textual.containers import Container, Vertical, ScrollableContainer
from textual import on
import json
from pathlib import Path
from rich.text import Text
import re

class AlignmentDetails(Static):
    """Displays details of the selected alignment node."""
    
    def parse_cigar(self, cigar: str) -> str:
        if not cigar: return "None"
        # Parse simple CIGAR like 100M1I50M
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
        if not ops: return cigar
        
        parsed = []
        for length, op in ops:
            desc = "Match" if op == 'M' else \
                   "Ins" if op == 'I' else \
                   "Del" if op == 'D' else \
                   "Skip" if op == 'N' else \
                   "Soft" if op == 'S' else \
                   "Hard" if op == 'H' else \
                   "Pad" if op == 'P' else \
                   "Equal" if op == '=' else \
                   "Mismatch" if op == 'X' else op
                   
            color = "green" if op in ['M', '='] else \
                    "red" if op in ['I', 'D', 'X'] else \
                    "yellow"
            
            parsed.append(f"[{color}]{desc}: {length}bp[/{color}]")
            
        return ", ".join(parsed)

    def update_details(self, data: dict):
        if not data:
            self.update("Select a node to view details.")
            return

        # Evidence nodes carry a "kind" discriminator so we can render the
        # right panel for each class of signal. Alignment nodes have "cigar".
        kind = data.get("_kind")
        if kind == "diploid":
            self._render_diploid(data)
            return
        if kind == "presence":
            self._render_presence(data)
            return
        if kind == "protein":
            self._render_protein(data)
            return
        if kind == "epistasis":
            self._render_epistasis(data)
            return
        if kind == "consequence":
            self._render_consequence(data)
            return

        if "cigar" not in data:
            self.update("Select an alignment or evidence node to view details.")
            return

        lines = []
        identity = data.get('identity', 0)
        color = "green" if identity >= 1.0 else "yellow" if identity > 0.95 else "red"
        lines.append(f"[bold cyan]Identity:[/bold cyan] [{color}]{identity:.4%}[/{color}]")
        lines.append(f"[bold cyan]Edit Distance (NM):[/bold cyan] {data.get('nm', '?')}")
        cigar = data.get('cigar', 'None')
        parsed_cigar = self.parse_cigar(cigar)
        lines.append(f"\n[bold cyan]CIGAR Summary:[/bold cyan]\n{parsed_cigar}")
        lines.append(f"\n[dim]Raw CIGAR: {cigar}[/dim]")
        self.update(Text.from_markup("\n".join(lines)))

    def _render_diploid(self, d: dict) -> None:
        lines = ["[bold cyan]Diploid Report[/bold cyan]\n"]
        for k in ("hap1_score", "hap2_score", "min_score", "max_score"):
            v = d.get(k)
            if isinstance(v, (int, float)):
                color = "green" if v >= 0.8 else "yellow" if v >= 0.5 else "red"
                lines.append(f"  [cyan]{k}:[/cyan] [{color}]{v:.3f}[/{color}]")
        if d.get("compound_het_lof"):
            lines.append("\n  [bold red]⚠ compound_het_lof = TRUE[/bold red]")
            lines.append("  [dim]both haplotypes independently loss-of-function[/dim]")
        constr = d.get("constraints") or {}
        if constr and any(v is not None for v in constr.values()):
            lines.append("\n[bold cyan]Gene constraint priors[/bold cyan]")
            for k in ("pli", "mis_z", "oe_lof", "clingen_haploinsufficiency"):
                if constr.get(k) is not None:
                    lines.append(f"  [cyan]{k}:[/cyan] {constr[k]}")
        self.update(Text.from_markup("\n".join(lines)))

    def _render_presence(self, d: dict) -> None:
        status = d.get("status", "?")
        color = {"intact": "green", "duplicated": "green",
                 "deleted": "red", "low_identity": "yellow",
                 "partial": "yellow", "uncertain": "yellow"}.get(status, "white")
        lines = [
            f"[bold cyan]Presence ({d.get('_hap', '?')})[/bold cyan]\n",
            f"  [cyan]status:[/cyan] [{color}]{status}[/{color}]",
            f"  [cyan]source:[/cyan] {d.get('source', '?')}",
        ]
        sid = d.get("sequence_identity")
        if isinstance(sid, (int, float)):
            lines.append(f"  [cyan]sequence_identity:[/cyan] {sid:.3f}")
        if d.get("copies", 1) and d.get("copies", 1) != 1:
            lines.append(f"  [cyan]copies:[/cyan] {d.get('copies')}")
        self.update(Text.from_markup("\n".join(lines)))

    def _render_protein(self, d: dict) -> None:
        ident = d.get("identity")
        color = "green" if (ident or 0) >= 0.99 else "yellow" if (ident or 0) >= 0.9 else "red"
        lines = [
            f"[bold cyan]Protein (T={d.get('transcript', '?')}, hap{d.get('haplotype', '?')})[/bold cyan]\n",
        ]
        if isinstance(ident, (int, float)):
            lines.append(f"  [cyan]identity:[/cyan] [{color}]{ident:.4f}[/{color}]")
        if d.get("ref_length") is not None:
            lines.append(f"  [cyan]ref_length:[/cyan] {d['ref_length']}")
        if d.get("hap_length") is not None:
            lines.append(f"  [cyan]hap_length:[/cyan] {d['hap_length']}")
        pts = d.get("premature_stop_at")
        if pts is not None:
            lines.append(f"  [bold red]premature stop at AA {pts}[/bold red]")
        fsr = d.get("frameshift_region") or {}
        if fsr:
            tag = "[green]restored[/green]" if fsr.get("restored") else "[red]not restored[/red]"
            lines.append(f"  [cyan]frameshift_region:[/cyan] {tag}")
        self.update(Text.from_markup("\n".join(lines)))

    def _render_epistasis(self, d: dict) -> None:
        resid = d.get("residual", 0)
        color = "red" if abs(resid) > 5 else "yellow" if abs(resid) > 2 else "white"
        flag = " [bold yellow][FLAGGED][/bold yellow]" if d.get("flagged") else ""
        lines = [
            f"[bold cyan]Epistasis residual (T={d.get('transcript', '?')}, hap{d.get('haplotype', '?')})[/bold cyan]\n",
            f"  [cyan]n_variants:[/cyan] {d.get('n_variants', '?')}",
            f"  [cyan]S_additive:[/cyan] {d.get('s_additive'):+.2f}",
            f"  [cyan]S_joint:[/cyan]    {d.get('s_joint'):+.2f}",
            f"  [cyan]residual:[/cyan]   [{color}]{resid:+.2f}[/{color}]{flag}",
            "",
            "[dim]positive residual → joint protein less deleterious than",
            "sum of per-variant log-odds (potential rescue)[/dim]",
        ]
        self.update(Text.from_markup("\n".join(lines)))

    def _render_consequence(self, d: dict) -> None:
        clnsig = d.get("clnsig") or ""
        clnsig_color = {
            "pathogenic": "red", "likely_pathogenic": "red",
            "benign": "green", "likely_benign": "green",
            "vus": "yellow", "conflicting": "yellow",
        }.get(clnsig, "white")
        lines = [
            f"[bold cyan]Consequence (hap{d.get('haplotype', '?')})[/bold cyan]\n",
            f"  [cyan]consequence:[/cyan] {d.get('consequence', '?')}",
        ]
        aa = d.get("amino_acid_change")
        if aa:
            lines.append(f"  [cyan]aa_change:[/cyan]   {aa}")
        if clnsig:
            lines.append(f"  [cyan]ClinVar:[/cyan]     [{clnsig_color}]{clnsig}[/{clnsig_color}]")
        self.update(Text.from_markup("\n".join(lines)))

class HapliExplorer(App):
    """A Textual app to explore Hapli alignment results."""

    CSS = """
    Screen {
        layout: vertical;
        background: $surface;
    }
    
    Header {
        dock: top;
        background: $accent;
        color: $text;
    }
    
    Footer {
        dock: bottom;
        background: $accent;
        color: $text;
    }
    
    #main-container {
        layout: horizontal;
        height: 1fr;
    }
    
    #left-pane {
        width: 40%;
        height: 100%;
        border-right: solid $primary;
        background: $surface-darken-1;
    }
    
    #right-pane {
        width: 60%;
        height: 100%;
        background: $surface;
        padding: 1 2;
    }
    
    Tree {
        background: $surface-darken-1;
        padding: 1;
        scrollbar-gutter: stable;
    }
    
    Tree:focus .tree--cursor {
        background: $primary;
        color: $text;
        text-style: bold;
    }
    
    Label {
        padding: 1;
        background: $primary-darken-2;
        color: $text;
        text-align: center;
        width: 100%;
    }
    """

    BINDINGS = [
        ("q", "quit", "Quit"),
        ("e", "expand_all", "Expand All"),
        ("c", "collapse_all", "Collapse All"),
    ]

    def __init__(self, json_path: Path):
        super().__init__()
        self.json_path = json_path
        self.data = self.load_data()

    def load_data(self):
        with open(self.json_path) as f:
            return json.load(f)

    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)
        with Container(id="main-container"):
            with Vertical(id="left-pane"):
                yield Label(f"Gene: {self.data.get('gene', 'Unknown')}")
                yield Tree("Transcripts")
            with Vertical(id="right-pane"):
                yield Label("Alignment Details")
                yield ScrollableContainer(AlignmentDetails(id="details"))
        yield Footer()

    def on_mount(self) -> None:
        tree = self.query_one(Tree)
        root = tree.root
        root.expand()
        ev = self.data.get("evidence") or {}

        # ── Evidence branch (schema v2) ────────────────────────────────────
        if ev:
            ev_node = root.add("[bold magenta]Evidence[/bold magenta]", expand=True)

            dip = ev.get("diploid") or {}
            if dip:
                s1 = dip.get("hap1_score"); s2 = dip.get("hap2_score")
                flag = " [bold red]⚠ COMPOUND HET LoF[/bold red]" if dip.get("compound_het_lof") else ""
                label = (
                    f"[bold]Diploid[/bold]: "
                    f"hap1={s1:.2f} hap2={s2:.2f}{flag}"
                    if isinstance(s1, (int, float)) and isinstance(s2, (int, float))
                    else "[bold]Diploid[/bold]"
                )
                dip_node = ev_node.add(label)
                dip_node.data = {**dip, "_kind": "diploid"}

            presence = ev.get("presence") or {}
            if presence:
                p_node = ev_node.add("[bold]Presence[/bold]", expand=True)
                for hap in ("hap1", "hap2"):
                    pc = presence.get(hap)
                    if isinstance(pc, dict):
                        status = pc.get("status", "?")
                        color = {"intact": "green", "duplicated": "green",
                                 "deleted": "red", "low_identity": "yellow",
                                 "partial": "yellow", "uncertain": "yellow"}.get(status, "white")
                        node = p_node.add(f"{hap}: [{color}]{status}[/{color}]")
                        node.data = {**pc, "_kind": "presence", "_hap": hap}

            prot = ev.get("protein") or []
            if prot:
                pr_node = ev_node.add(f"[bold]Protein diffs[/bold] ({len(prot)})", expand=True)
                for p in prot:
                    ident = p.get("identity", 0) or 0
                    color = "green" if ident >= 0.99 else "yellow" if ident >= 0.9 else "red"
                    tag = ""
                    if p.get("premature_stop_at") is not None:
                        tag = " [red]PTC[/red]"
                    elif (p.get("frameshift_region") or {}).get("restored"):
                        tag = " [green]frameshift restored[/green]"
                    node = pr_node.add(
                        f"T={p.get('transcript', '?')} hap{p.get('haplotype', '?')}: "
                        f"[{color}]{ident:.3f}[/{color}]{tag}"
                    )
                    node.data = {**p, "_kind": "protein"}

            epi = ev.get("epistasis") or []
            epi_with_residual = [r for r in epi if r.get("residual") is not None]
            if epi_with_residual:
                ep_node = ev_node.add(f"[bold]Epistasis[/bold] ({len(epi_with_residual)})", expand=True)
                for r in epi_with_residual:
                    resid = r.get("residual", 0)
                    color = "red" if abs(resid) > 5 else "yellow" if abs(resid) > 2 else "white"
                    flag = " [bold yellow]⚑[/bold yellow]" if r.get("flagged") else ""
                    node = ep_node.add(
                        f"T={r.get('transcript', '?')} hap{r.get('haplotype', '?')}: "
                        f"residual=[{color}]{resid:+.2f}[/{color}]{flag}"
                    )
                    node.data = {**r, "_kind": "epistasis"}

            cons = ev.get("consequence") or []
            if cons:
                c_node = ev_node.add(f"[bold]Consequences[/bold] ({len(cons)})")
                for c in cons:
                    clnsig = c.get("clnsig") or ""
                    tag = ""
                    if clnsig:
                        cs_color = {"pathogenic": "red", "likely_pathogenic": "red",
                                    "benign": "green", "likely_benign": "green",
                                    "vus": "yellow"}.get(clnsig, "white")
                        tag = f" [{cs_color}]({clnsig})[/{cs_color}]"
                    node = c_node.add(
                        f"hap{c.get('haplotype', '?')} {c.get('consequence', '?')}{tag}"
                    )
                    node.data = {**c, "_kind": "consequence"}

        # ── Transcripts branch (alignment view, schema v1 compatible) ──────
        transcripts = self.data.get("transcripts", [])
        if transcripts:
            tx_node = root.add("[bold]Transcripts[/bold]", expand=True)
            for t in transcripts:
                t_node = tx_node.add(f"[bold]{t['id']}[/bold]", expand=True)
                alignments = t.get("alignments", {})
                for hap, info in alignments.items():
                    if isinstance(info, str):
                        t_node.add(f"{hap}: [red]{info}[/red]")
                        continue
                    identity = info.get('identity', 0)
                    color = "green" if identity >= 1.0 else "yellow" if identity > 0.95 else "red"
                    hap_node = t_node.add(f"{hap}: [{color}]{identity:.2%}[/{color}]")
                    hap_node.data = info

    @on(Tree.NodeSelected)
    def show_details(self, event: Tree.NodeSelected) -> None:
        node = event.node
        data = node.data
        if data:
            self.query_one("#details", AlignmentDetails).update_details(data)

    def action_expand_all(self):
        tree = self.query_one(Tree)
        tree.root.expand_all()

    def action_collapse_all(self):
        tree = self.query_one(Tree)
        tree.root.collapse_all()
        tree.root.expand() # Keep root open
