from textual.app import App, ComposeResult
from textual.widgets import Header, Footer, Tree, Static, Label, ScrollableContainer
from textual.containers import Container, Horizontal, Vertical
from textual.widgets.tree import TreeNode
from textual import on
import json
from pathlib import Path
from rich.text import Text
from rich.syntax import Syntax
from rich.panel import Panel
from rich.layout import Layout
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
            self.update("Select a feature to view details.")
            return

        # Header
        lines = []
        lines.append(f"[bold cyan]Feature:[/bold cyan] {data.get('feature_type', 'unknown')} {data.get('feature_id', 'unknown')}")
        lines.append(f"[bold cyan]Region:[/bold cyan] {data.get('target_start', '?')}-{data.get('target_end', '?')}")
        
        identity = data.get('identity', 0)
        color = "green" if identity >= 1.0 else "yellow" if identity > 0.9 else "red"
        lines.append(f"[bold cyan]Identity:[/bold cyan] [{color}]{identity:.4%}[/{color}]")
        lines.append(f"[bold cyan]MAPQ:[/bold cyan] {data.get('mapq', 0)})
        
        cigar = data.get('cigar', 'None')
        parsed_cigar = self.parse_cigar(cigar)
        lines.append(f"\n[bold cyan]CIGAR Summary:[/bold cyan]\n{parsed_cigar}")
        lines.append(f"\n[dim]Raw CIGAR: {cigar}[/dim]")
        
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

    .node-perfect {
        color: $success;
    }
    
    .node-good {
        color: $warning;
    }
    
    .node-poor {
        color: $error;
    }
    
    Label {
        padding: 1;
        background: $primary-darken-2;
        color: $text;
        text-align: center;
        width: 100%;
    }
    
    #details {
        padding-top: 1;
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
                yield Label("Alignment Hierarchy")
                yield Tree("Haplotypes")
            with Vertical(id="right-pane"):
                yield Label("Feature Details")
                yield ScrollableContainer(AlignmentDetails(id="details"))
        yield Footer()

    def on_mount(self) -> None:
        tree = self.query_one(Tree)
        root = tree.root
        root.expand()
        
        for hap_name, hap_data in self.data.items():
            # Haplotype Node
            hap_node = root.add(f"[bold]{hap_name}[/bold]", expand=True)
            hap_node.data = hap_data
            
            # Recursively add children
            self.add_children(hap_node, hap_data)

    def add_children(self, parent_node: TreeNode, data: dict):
        children = data.get("children", [])
        # Sort children by start position for better readability
        children.sort(key=lambda x: x.get('target_start', 0))
        
        for child in children:
            feat_type = child['feature_type']
            feat_id = child['feature_id']
            identity = child.get('identity', 0)
            
            # Icon based on type
            icon = "🧬" if feat_type == "gene" else \
                   "📄" if feat_type == "transcript" else \
                   "🧱" if feat_type == "exon" else \
                   "⚙️" if feat_type == "CDS" else "•"

            label = f"{icon} {feat_type}: {feat_id}"
            
            node = parent_node.add(label, expand=False)
            node.data = child
            
            # Style the label based on identity (Textual supports some rich styling in tree labels, 
            # but class-based styling is more robust for rows. 
            # However, TreeNode doesn't easily accept CSS classes for just the label part unless customized.
            # We can use Rich markup in the label.)
            
            if identity >= 1.0:
                node.label = f"[green]{label}[/green]"
            elif identity > 0.95:
                node.label = f"[yellow]{label}[/yellow]"
            else:
                node.label = f"[red]{label}[/red]"
            
            self.add_children(node, child)

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

if __name__ == "__main__":
    import sys
    app = HapliExplorer(Path(sys.argv[1]))
    app.run()