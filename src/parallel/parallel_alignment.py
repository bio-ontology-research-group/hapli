"""
Parallel alignment module.

Provides functionality to parallelize minimap2 executions across multiple 
processes with smart batching and result aggregation.
"""

import logging
from typing import List, Dict, Any, Tuple, Optional, Union, Set
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from src.parallel.task_manager import execute_parallel, create_worker_pool, TaskChunker
from src.alignment.minimap_wrapper import MinimapAligner

logger = logging.getLogger(__name__)

class BatchAlignmentTask:
    """
    Represents a batch of sequences to be aligned.
    
    This class encapsulates all data needed for a worker process
    to perform alignment tasks independently.
    """
    
    def __init__(self, 
                features: List[SeqFeature],
                reference_seq: SeqRecord,
                aligner_params: Dict[str, Any] = None):
        """
        Initialize a batch alignment task.
        
        Args:
            features: List of features to align
            reference_seq: Reference sequence to align against
            aligner_params: Parameters for the minimap aligner
        """
        self.features = features
        self.reference_seq = reference_seq
        self.aligner_params = aligner_params or {}
        
    def __len__(self) -> int:
        """Get the number of features in this batch."""
        return len(self.features)


def _process_alignment_batch(batch: BatchAlignmentTask) -> List[Tuple[str, List]]:
    """
    Process a batch of alignments.
    
    Args:
        batch: BatchAlignmentTask containing features and reference sequence
        
    Returns:
        List of tuples (feature_id, alignment_results)
    """
    results = []
    
    # Create a new aligner instance for this batch
    aligner = MinimapAligner(**batch.aligner_params)
    aligner.load_reference(batch.reference_seq)
    
    # Process each feature in the batch
    for feature in batch.features:
        feature_id = feature.id if hasattr(feature, 'id') else str(feature.qualifiers.get('ID', [None])[0])
        if not feature_id:
            logger.warning(f"Feature without ID encountered: {feature}")
            continue
            
        try:
            # Align the feature
            alignments = aligner.align_feature(feature, batch.reference_seq)
            results.append((feature_id, alignments))
        except Exception as e:
            logger.error(f"Error aligning feature {feature_id}: {str(e)}")
            results.append((feature_id, []))
            
    return results


class ParallelAligner:
    """
    Manages parallel alignment of features using minimap2.
    
    This class handles:
    - Distributing alignment tasks across multiple processes
    - Smart batching to reduce overhead
    - Result aggregation
    - Progress tracking
    """
    
    def __init__(self, 
                num_workers: Optional[int] = None,
                batch_size: int = 50,
                pool_type: str = 'process',
                minimap_preset: str = 'splice',
                **minimap_params):
        """
        Initialize parallel aligner.
        
        Args:
            num_workers: Number of worker processes
            batch_size: Number of features per batch
            pool_type: Type of worker pool ('process' or 'thread')
            minimap_preset: Preset for minimap2 
            **minimap_params: Additional parameters for minimap2
        """
        self.num_workers = num_workers
        self.batch_size = batch_size
        self.pool_type = pool_type
        self.minimap_params = {'preset': minimap_preset, **minimap_params}
        
    def align_features(self, 
                      features: List[SeqFeature],
                      reference_seq: SeqRecord) -> Dict[str, List]:
        """
        Align multiple features in parallel.
        
        Args:
            features: List of features to align
            reference_seq: Reference sequence to align against
            
        Returns:
            Dictionary mapping feature IDs to alignment results
        """
        if not features:
            logger.warning("No features provided for alignment")
            return {}
            
        logger.info(f"Aligning {len(features)} features using {self.num_workers} workers")
        
        # Create batches
        batches = self._create_batches(features, reference_seq)
        logger.info(f"Created {len(batches)} batches (batch size: {self.batch_size})")
        
        # Execute alignment in parallel
        with create_worker_pool(self.pool_type, self.num_workers) as pool:
            batch_results = pool.map(_process_alignment_batch, batches)
            
        # Aggregate results
        all_results = {}
        for batch_result in batch_results:
            for feature_id, alignments in batch_result:
                all_results[feature_id] = alignments
                
        logger.info(f"Completed alignment of {len(all_results)} features")
        return all_results
        
    def _create_batches(self, 
                       features: List[SeqFeature], 
                       reference_seq: SeqRecord) -> List[BatchAlignmentTask]:
        """
        Create batches of features for parallel processing.
        
        Args:
            features: List of features to align
            reference_seq: Reference sequence
            
        Returns:
            List of BatchAlignmentTask objects
        """
        # Group features into batches
        feature_batches = TaskChunker.chunk_by_size(features, self.batch_size)
        
        # Create task objects for each batch
        return [
            BatchAlignmentTask(
                features=batch,
                reference_seq=reference_seq,
                aligner_params=self.minimap_params
            )
            for batch in feature_batches
        ]
