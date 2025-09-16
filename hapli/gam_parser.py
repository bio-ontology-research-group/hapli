#!/usr/bin/env python3
"""
GAM file parser for extracting alignment data and organizing by sample/haplotype.

This module provides functionality to parse GAM (Graph Alignment/Map) files
and organize the alignment data by sample and haplotype path names.
"""

import json
import logging
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import re

logger = logging.getLogger(__name__)


class GAMParser:
    """
    Parse GAM files and organize alignments by sample/haplotype path names.
    
    This class reads GAM alignments using vg tools, extracts sample and haplotype
    information from path names, and builds a hierarchical feature structure.
    """
    
    def __init__(self, gam_file: Path, vg_executable: str = "vg"):
        """
        Initialize GAM parser.
        
        Args:
            gam_file: Path to the GAM file to parse
            vg_executable: Path to vg executable (default: "vg")
        """
        self.gam_file = Path(gam_file)
        self.vg_executable = vg_executable
        self.alignments = []
        self._validate_inputs()
    
    def _validate_inputs(self) -> None:
        """Validate input files and executables."""
        if not self.gam_file.exists():
            raise FileNotFoundError(f"GAM file not found: {self.gam_file}")
        
        # Check if vg is available
        try:
            result = subprocess.run(
                [self.vg_executable, "version"], 
                capture_output=True, check=True, text=True
            )
            logger.debug(f"vg version check successful: {result.stdout.strip()}")
        except FileNotFoundError:
            raise RuntimeError(f"vg executable not found at the specified path: {self.vg_executable}")
        except subprocess.CalledProcessError as e:
            error_message = (
                f"vg executable at '{self.vg_executable}' failed to run. "
                f"Return code: {e.returncode}\n"
                f"Stderr: {e.stderr.strip()}"
            )
            raise RuntimeError(error_message)
    
    def _run_vg_command(self, command: List[str]) -> str:
        """
        Run a vg command and return output.
        
        Args:
            command: List of command arguments
            
        Returns:
            Command output as string
        """
        full_command = [self.vg_executable] + command
        try:
            result = subprocess.run(full_command, capture_output=True, 
                                  text=True, check=True)
            return result.stdout
        except subprocess.CalledProcessError as e:
            logger.error(f"vg command failed: {' '.join(full_command)}")
            logger.error(f"Error: {e.stderr}")
            raise
    
    def load_alignments(self) -> None:
        """Load GAM alignments using vg view to convert to JSON."""
        logger.info(f"Loading alignments from {self.gam_file}")
        
        # Convert GAM to JSON using vg view
        json_output = self._run_vg_command(["view", "-a", str(self.gam_file)])
        
        # Parse JSON lines
        self.alignments = []
        for line in json_output.strip().split('\n'):
            if line.strip():
                try:
                    alignment = json.loads(line)
                    self.alignments.append(alignment)
                except json.JSONDecodeError as e:
                    logger.warning(f"Failed to parse alignment JSON: {e}")
                    continue
        
        logger.info(f"Loaded {len(self.alignments)} alignments")
    
    def _extract_sample_haplotype_from_path(self, path_name: str) -> Tuple[str, str]:
        """
        Extract sample and haplotype information from path name.
        
        Common patterns:
        - sample_name#0#chr1 -> sample: sample_name, haplotype: 0
        - sample_name.1.chr1 -> sample: sample_name, haplotype: 1
        - GRCh38#0#chr1 -> sample: GRCh38, haplotype: 0
        
        Args:
            path_name: Path name from alignment
            
        Returns:
            Tuple of (sample_name, haplotype)
        """
        # Pattern 1: sample#haplotype#chromosome
        pattern1 = re.match(r'^([^#]+)#(\d+)#', path_name)
        if pattern1:
            return pattern1.group(1), pattern1.group(2)
        
        # Pattern 2: sample.haplotype.chromosome
        pattern2 = re.match(r'^([^.]+)\.(\d+)\.', path_name)
        if pattern2:
            return pattern2.group(1), pattern2.group(2)
        
        # Pattern 3: sample_haplotype_chromosome (underscore separated)
        pattern3 = re.match(r'^(.+)_(\d+)_', path_name)
        if pattern3:
            return pattern3.group(1), pattern3.group(2)
        
        # Default: treat entire path as sample, haplotype as "unknown"
        return path_name, "unknown"
    
    def _extract_feature_info_from_alignment(self, alignment: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract feature information from alignment metadata.
        
        Args:
            alignment: Alignment dictionary from GAM
            
        Returns:
            Dictionary containing feature information
        """
        feature_info = {
            'read_name': alignment.get('name', 'unknown'),
            'sequence': alignment.get('sequence', ''),
            'quality': alignment.get('quality', ''),
            'mapping_quality': None,
            'is_reverse': False,
            'path_positions': [],
            'score': alignment.get('score', 0),
            'identity': alignment.get('identity', 0.0)
        }
        
        # Extract mapping information
        if 'path' in alignment:
            path_info = alignment['path']
            feature_info['mapping_quality'] = path_info.get('mapping_quality', 0)
            
            # Extract path positions from mappings
            if 'mapping' in path_info:
                for mapping in path_info['mapping']:
                    position_info = {
                        'node_id': mapping.get('position', {}).get('node_id', 0),
                        'offset': mapping.get('position', {}).get('offset', 0),
                        'is_reverse': mapping.get('position', {}).get('is_reverse', False)
                    }
                    feature_info['path_positions'].append(position_info)
                    
                    # Update reverse flag if any mapping is reverse
                    if position_info['is_reverse']:
                        feature_info['is_reverse'] = True
        
        # Try to extract feature type from read name
        feature_info['feature_type'] = self._infer_feature_type(feature_info['read_name'])
        
        return feature_info
    
    def _infer_feature_type(self, read_name: str) -> str:
        """
        Infer feature type from read name patterns.
        
        Args:
            read_name: Name of the read/feature
            
        Returns:
            Inferred feature type
        """
        read_name_lower = read_name.lower()
        
        if 'gene' in read_name_lower:
            return 'gene'
        elif 'exon' in read_name_lower:
            return 'exon'
        elif 'cds' in read_name_lower:
            return 'CDS'
        elif 'utr' in read_name_lower:
            return 'UTR'
        elif 'transcript' in read_name_lower or 'mrna' in read_name_lower:
            return 'mRNA'
        else:
            return 'feature'
    
    def _build_feature_hierarchy(self, features: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Build hierarchical structure from flat list of features.
        
        Args:
            features: List of feature dictionaries
            
        Returns:
            Hierarchical feature structure
        """
        hierarchy = defaultdict(list)
        
        # Group by feature type
        for feature in features:
            feature_type = feature['feature_type']
            hierarchy[feature_type].append(feature)
        
        # Sort features within each type by score (descending)
        for feature_type in hierarchy:
            hierarchy[feature_type].sort(key=lambda x: x['score'], reverse=True)
        
        return dict(hierarchy)
    
    def group_alignments_by_sample_haplotype(self) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """
        Group alignments by sample and haplotype path names.
        
        Returns:
            Nested dictionary: {sample: {haplotype: {feature_type: [alignment_data]}}}
        """
        if not self.alignments:
            self.load_alignments()
        
        grouped_data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        
        for alignment in self.alignments:
            # Get path name from alignment
            path_name = None
            if 'path' in alignment and 'name' in alignment['path']:
                path_name = alignment['path']['name']
            elif 'refpos' in alignment:
                # Alternative path name location
                for refpos in alignment['refpos']:
                    if 'name' in refpos:
                        path_name = refpos['name']
                        break
            
            if not path_name:
                logger.warning("No path name found in alignment, skipping")
                continue
            
            # Extract sample and haplotype
            sample, haplotype = self._extract_sample_haplotype_from_path(path_name)
            
            # Extract feature information
            feature_info = self._extract_feature_info_from_alignment(alignment)
            feature_type = feature_info['feature_type']
            
            # Add to grouped data
            grouped_data[sample][haplotype][feature_type].append(feature_info)
        
        # Build hierarchies for each sample/haplotype combination
        result = {}
        for sample in grouped_data:
            result[sample] = {}
            for haplotype in grouped_data[sample]:
                # Flatten all features for this sample/haplotype
                all_features = []
                for feature_type, features in grouped_data[sample][haplotype].items():
                    all_features.extend(features)
                
                # Build hierarchy
                result[sample][haplotype] = self._build_feature_hierarchy(all_features)
        
        return result
    
    def get_alignment_statistics(self) -> Dict[str, Any]:
        """
        Get statistics about the loaded alignments.
        
        Returns:
            Dictionary containing alignment statistics
        """
        if not self.alignments:
            self.load_alignments()
        
        stats = {
            'total_alignments': len(self.alignments),
            'samples': set(),
            'haplotypes': set(),
            'feature_types': defaultdict(int),
            'avg_score': 0.0,
            'avg_identity': 0.0
        }
        
        total_score = 0
        total_identity = 0
        
        for alignment in self.alignments:
            # Extract path info for sample/haplotype counting
            path_name = None
            if 'path' in alignment and 'name' in alignment['path']:
                path_name = alignment['path']['name']
            
            if path_name:
                sample, haplotype = self._extract_sample_haplotype_from_path(path_name)
                stats['samples'].add(sample)
                stats['haplotypes'].add(haplotype)
            
            # Feature type counting
            feature_info = self._extract_feature_info_from_alignment(alignment)
            stats['feature_types'][feature_info['feature_type']] += 1
            
            # Score and identity accumulation
            total_score += alignment.get('score', 0)
            total_identity += alignment.get('identity', 0.0)
        
        # Calculate averages
        if stats['total_alignments'] > 0:
            stats['avg_score'] = total_score / stats['total_alignments']
            stats['avg_identity'] = total_identity / stats['total_alignments']
        
        # Convert sets to lists for JSON serialization
        stats['samples'] = list(stats['samples'])
        stats['haplotypes'] = list(stats['haplotypes'])
        stats['feature_types'] = dict(stats['feature_types'])
        
        return stats
    
    def save_grouped_data(self, output_file: Path, 
                         grouped_data: Optional[Dict[str, Dict[str, Dict[str, Any]]]] = None) -> None:
        """
        Save grouped alignment data to JSON file.
        
        Args:
            output_file: Output JSON file path
            grouped_data: Grouped data to save (if None, will group alignments first)
        """
        if grouped_data is None:
            grouped_data = self.group_alignments_by_sample_haplotype()
        
        with open(output_file, 'w') as f:
            json.dump(grouped_data, f, indent=2)
        
        logger.info(f"Saved grouped alignment data to {output_file}")


def main():
    """Command-line interface for GAM parser."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Parse GAM files and group by sample/haplotype")
    parser.add_argument("gam_file", help="Input GAM file")
    parser.add_argument("-o", "--output", help="Output JSON file for grouped data")
    parser.add_argument("--stats", help="Output file for alignment statistics")
    parser.add_argument("--vg", default="vg", help="Path to vg executable")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Parse GAM file
    parser_obj = GAMParser(Path(args.gam_file), args.vg)
    
    # Group alignments
    grouped_data = parser_obj.group_alignments_by_sample_haplotype()
    
    # Save grouped data
    if args.output:
        parser_obj.save_grouped_data(Path(args.output), grouped_data)
    
    # Save statistics
    if args.stats:
        stats = parser_obj.get_alignment_statistics()
        with open(args.stats, 'w') as f:
            json.dump(stats, f, indent=2)
    
    # Print summary
    stats = parser_obj.get_alignment_statistics()
    print(f"Processed {stats['total_alignments']} alignments")
    print(f"Found {len(stats['samples'])} samples and {len(stats['haplotypes'])} haplotypes")
    print(f"Feature types: {list(stats['feature_types'].keys())}")


if __name__ == "__main__":
    main()
