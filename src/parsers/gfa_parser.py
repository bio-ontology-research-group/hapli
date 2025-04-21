"""
Parser for Graphical Fragment Assembly (GFA) files.
Uses the PyGFA library to handle GFA format parsing.
"""
import os
import logging
from typing import Dict, List, Optional, Set, Tuple, Any
import pygfa
from pygfa.graph_element.parser import header, segment, link, path, line

class GFAParser:
    """
    Parser for Graphical Fragment Assembly (GFA) files.
    Handles both GFA1 and GFA2 formats through PyGFA.
    """
    
    def __init__(self):
        """Initialize the GFA parser."""
        self.logger = logging.getLogger(__name__)
        self.gfa = None
        self._segments = {}
        self._paths = {}

    def parse(self, filepath: str) -> pygfa.graph.GraphicalFragmentAssembly:
        """
        Parse a GFA file and return the PyGFA graph object.
        
        Args:
            filepath: Path to the GFA file
            
        Returns:
            A PyGFA GraphicalFragmentAssembly object
            
        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file is empty or malformed
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"GFA file not found: {filepath}")
            
        self.logger.info(f"Parsing GFA file: {filepath}")
        
        try:
            # Create a new GFA object
            self.gfa = pygfa.graph.GraphicalFragmentAssembly()
            
            # Parse the file line by line to handle potential format issues
            with open(filepath, 'r') as gfa_file:
                line_count = 0
                for line_num, line_str in enumerate(gfa_file, 1):
                    line_str = line_str.strip()
                    if not line_str or line_str.startswith('#'):
                        continue
                        
                    line_count += 1
                    fields = line_str.split('\t')
                    if not fields:
                        continue
                        
                    # Process based on record type
                    record_type = fields[0]
                    try:
                        if record_type == 'H':
                            header_line = header.Header.from_string(line_str)
                            self.gfa.add_header(header_line)
                        elif record_type == 'S':
                            seg = segment.SegmentV1.from_string(line_str)
                            self.gfa.add_segment(seg)
                            self._segments[seg.name] = seg
                        elif record_type == 'L':
                            link_line = link.Link.from_string(line_str)
                            self.gfa.add_edge(link_line)
                        elif record_type == 'P':
                            path_line = path.Path.from_string(line_str)
                            self.gfa.add_path(path_line)
                            self._paths[path_line.name] = path_line
                        # Add support for other record types as needed
                    except Exception as e:
                        self.logger.warning(f"Error parsing line {line_num}: {e}")
                
            if line_count == 0:
                raise ValueError("Empty or invalid GFA file")
                
            self.logger.info(f"Successfully parsed GFA with {len(self._segments)} segments and {len(self._paths)} paths")
            return self.gfa
            
        except Exception as e:
            self.logger.error(f"Failed to parse GFA file: {e}")
            raise
            
    def get_segments(self) -> Dict[str, Any]:
        """Return all segments in the GFA file."""
        return self._segments
        
    def get_paths(self) -> Dict[str, Any]:
        """Return all paths in the GFA file."""
        return self._paths
    
    def get_segment_sequence(self, segment_id: str) -> Optional[str]:
        """Get the sequence for a specific segment."""
        if not self.gfa or segment_id not in self._segments:
            return None
        return self._segments[segment_id].sequence
