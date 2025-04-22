"""
Parallel processing module.

Provides functionality for parallel execution using processes and threads,
with support for hierarchical dependencies and efficient task chunking.
"""

from src.parallel.task_manager import (
    create_worker_pool,
    execute_parallel,
    TaskChunker,
    ProgressTracker,
    ProcessWorkerPool,
    ThreadWorkerPool
)

from src.parallel.parallel_alignment import (
    ParallelAligner,
    BatchAlignmentTask
)

from src.parallel.hierarchical_executor import (
    HierarchicalExecutor,
    Task,
    execute_hierarchical_tasks
)
