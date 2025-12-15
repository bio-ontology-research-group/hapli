from textual.app import App, ComposeResult
from textual.widgets import Header, Footer, Tree, Static, Label
from textual.containers import Container, Horizontal, Vertical
from textual.widgets.tree import TreeNode
from textual import on
import json
from pathlib import Path
from rich.text import Text
from rich.syntax import Syntax

class AlignmentDetails(Static):
    """Displays details of the selected alignment node."""
    
    def update_details(self, data: dict):
        if not data:
            self.update("Select a feature to view details.")
            return

        lines = []
        lines.append(f"[b]Feature:[/b] {data.get('feature_type', 'unknown')} {data.get('feature_id', 'unknown')}")
        lines.append(f"[b]Region:[/b] {data.get('target_start', '?')}-{data.get('target_end', '?')}")
        
        identity = data.get('identity', 0)
        color = "green" if identity >= 1.0 else "yellow" if identity > 0.9 else "red"
        lines.append(f"[b]Identity:[/b] [{color}]{identity:.4f}[/{color}]")
        
        lines.append(f"[b]MAPQ:[/b] {data.get('mapq', 0)}")
        lines.append(f"[b]CIGAR:[/b] {data.get('cigar', 'None')}")
        
        if 'children' in data and data['children']:
            lines.append(f"\n[b]Children:[/b] {len(data['children'])}")
        
        self.update(Text.from_markup("\n".join(lines)))

class HapliExplorer(App):
    """A Textual app to explore Hapli alignment results."""

    CSS = """
    Screen {
        layout: vertical;
    }
    
    Header {
        dock: top;
    }
    
    Footer {
        dock: bottom;
    }
    
    #main-container {
        layout: horizontal;
        height: 1fr;
    }
    
    #tree-container {
        width: 40%;
        height: 100%;
        border: solid green;
    }
    
    #details-container {
        width: 60%;
        height: 100%;
        border: solid blue;
        padding: 1 2;
    }
    
    Tree {
        padding: 1;
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
        yield Header()
        with Container(id="main-container"):
            with Container(id="tree-container"):
                yield Label("[b]Alignment Hierarchy[/b]")
                yield Tree("Haplotypes")
            with Vertical(id="details-container"):
                yield Label("[b]Feature Details[/b]")
                yield AlignmentDetails(id="details")
        yield Footer()

    def on_mount(self) -> None:
        tree = self.query_one(Tree)
        root = tree.root
        root.expand()
        
        for hap_name, hap_data in self.data.items():
            # Haplotype Node
            hap_node = root.add(f"[b]{hap_name}[/b]", expand=True)
            hap_node.data = hap_data
            
            # Recursively add children
            self.add_children(hap_node, hap_data)

    def add_children(self, parent_node: TreeNode, data: dict):
        # The top level dict from JSON is just the root object, but its children are in 'children' list
        # If this is the haplotype root, it IS the gene result object usually.
        
        children = data.get("children", [])
        for child in children:
            label = f"{child['feature_type']}: {child['feature_id']}"
            
            # Color code by identity
            identity = child.get('identity', 0)
            if identity < 1.0:
                label = f"[yellow]{label}[/yellow]"
            
            node = parent_node.add(label, expand=False)
            node.data = child
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
        tree.root.expand()

if __name__ == "__main__":
    import sys
    app = HapliExplorer(Path(sys.argv[1]))
    app.run()
