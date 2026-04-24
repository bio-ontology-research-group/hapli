import subprocess
import tempfile
from pathlib import Path
import pysam
import logging
from typing import List, Optional

class SequenceAligner:
    def __init__(self, minimap2_path: str = "minimap2"):
        self.minimap2_path = minimap2_path
        self.logger = logging.getLogger(__name__)

    def align(self, query_seq: str, query_name: str, target_seq: str, target_name: str) -> List[pysam.AlignedSegment]:
        """
        Aligns query sequence to target sequence using minimap2.
        Returns list of pysam.AlignedSegment.
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix=".fa", delete=False) as target_f, \
             tempfile.NamedTemporaryFile(mode='w', suffix=".fa", delete=False) as query_f:
            
            target_f.write(f">{target_name}\n{target_seq}\n")
            target_path = target_f.name
            
            query_f.write(f">{query_name}\n{query_seq}\n")
            query_path = query_f.name

        try:
            # Run minimap2
            # -a: output SAM
            # -x splice: splice-aware alignment (good for transcripts mapping to genome)
            # --cs: output cs tag for precise difference (variant calling)
            cmd = [self.minimap2_path, "-a", "-x", "splice", "--cs", target_path, query_path]
            
            self.logger.debug(f"Running command: {' '.join(cmd)}")
            
            # Write output directly to a temp SAM file
            with tempfile.NamedTemporaryFile(mode='w', suffix=".sam", delete=False) as sam_out:
                subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, check=True, text=True)
                sam_path = sam_out.name
            
            # Parse SAM output
            alignments = []
            try:
                with pysam.AlignmentFile(sam_path, "r", check_sq=False) as sam_file:
                    for record in sam_file:
                        if not record.is_unmapped:
                            alignments.append(record)
            finally:
                Path(sam_path).unlink(missing_ok=True)
            
            alignments.sort(key=lambda x: x.mapping_quality, reverse=True)
            return alignments

        except subprocess.CalledProcessError as e:
            self.logger.error(f"minimap2 failed: {e.stderr}")
            raise
        except FileNotFoundError:
             self.logger.error("minimap2 not found. Please ensure it is installed and in your PATH.")
             raise
        finally:
            Path(target_path).unlink(missing_ok=True)
            Path(query_path).unlink(missing_ok=True)
