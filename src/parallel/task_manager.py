"""
Task manager module for parallel processing.

Provides worker pool implementations using both multiprocessing and threading.
Includes task chunking, progress tracking, and proper resource management.
"""

import os
import time
import logging
import multiprocessing
from multiprocessing import Pool as ProcessPool
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import List, Callable, Any, Dict, Union, Optional, Tuple, Generator, Iterable

logger = logging.getLogger(__name__)

class TaskChunker:
    """Utility class to divide a list of tasks into chunks for parallel processing."""
    
    @staticmethod
    def chunk_tasks(tasks: List[Any], num_chunks: int) -> List[List[Any]]:
        """
        Divide tasks into a specified number of chunks.
        
        Args:
            tasks: List of tasks to divide
            num_chunks: Number of chunks to divide into
            
        Returns:
            List of task chunks
        """
        if not tasks:
            return []
            
        # Adjust num_chunks if it's greater than the number of tasks
        num_chunks = min(num_chunks, len(tasks))
        
        # Calculate approximate size of each chunk
        chunk_size = len(tasks) // num_chunks
        remainder = len(tasks) % num_chunks
        
        chunks = []
        start = 0
        for i in range(num_chunks):
            # Add one extra item to the first 'remainder' chunks
            end = start + chunk_size + (1 if i < remainder else 0)
            chunks.append(tasks[start:end])
            start = end
            
        return chunks
    
    @staticmethod
    def chunk_by_size(tasks: List[Any], chunk_size: int) -> List[List[Any]]:
        """
        Divide tasks into chunks of specified size.
        
        Args:
            tasks: List of tasks to divide
            chunk_size: Maximum size of each chunk
            
        Returns:
            List of task chunks
        """
        if not tasks:
            return []
            
        return [tasks[i:i + chunk_size] for i in range(0, len(tasks), chunk_size)]


class ProgressTracker:
    """Tracks and reports progress of parallel tasks."""
    
    def __init__(self, total_tasks: int, report_interval: float = 2.0):
        """
        Initialize progress tracker.
        
        Args:
            total_tasks: Total number of tasks to track
            report_interval: Minimum interval between progress reports in seconds
        """
        self.total_tasks = total_tasks
        
        # Use Value from multiprocessing for thread-safe counter that works across processes
        try:
            self.completed_tasks = multiprocessing.Value('i', 0)
            self.is_mp_value = True
        except:
            # Fallback if multiprocessing Value fails
            self.completed_tasks = 0
            self.is_mp_value = False
            
        self.start_time = time.time()
        self.report_interval = report_interval
        self.last_report_time = self.start_time
        self.lock = threading.RLock()  # Use RLock for thread safety
    
    def increment(self, count: int = 1) -> None:
        """
        Increment the number of completed tasks.
        
        Args:
            count: Number of tasks completed
        """
        with self.lock:
            if self.is_mp_value:
                with self.completed_tasks.get_lock():
                    self.completed_tasks.value += count
                completed = self.completed_tasks.value
            else:
                self.completed_tasks += count
                completed = self.completed_tasks
            
            current_time = time.time()
            if (current_time - self.last_report_time >= self.report_interval or 
                completed >= self.total_tasks):
                self._report_progress()
                self.last_report_time = current_time
    
    def _report_progress(self) -> None:
        """Report current progress."""
        if self.is_mp_value:
            with self.completed_tasks.get_lock():
                completed = self.completed_tasks.value
        else:
            completed = self.completed_tasks
            
        percent = (completed / self.total_tasks) * 100 if self.total_tasks > 0 else 0
        elapsed = time.time() - self.start_time
        
        # Calculate estimated time remaining
        if completed > 0:
            tasks_per_second = completed / elapsed
            remaining_tasks = self.total_tasks - completed
            eta = remaining_tasks / tasks_per_second if tasks_per_second > 0 else 0
            eta_str = f", ETA: {eta:.1f}s" if eta > 0 else ""
        else:
            eta_str = ""
            
        logger.info(f"Progress: {completed}/{self.total_tasks} ({percent:.1f}%) "
                   f"in {elapsed:.1f}s{eta_str}")
    
    def get_completed(self) -> int:
        """
        Get the number of completed tasks.
        
        Returns:
            Number of completed tasks
        """
        if self.is_mp_value:
            with self.completed_tasks.get_lock():
                return self.completed_tasks.value
        else:
            return self.completed_tasks


class BaseWorkerPool:
    """Base class for worker pools."""
    
    def __init__(self, num_workers: Optional[int] = None, 
                 track_progress: bool = True):
        """
        Initialize worker pool.
        
        Args:
            num_workers: Number of worker processes/threads
            track_progress: Whether to track and report progress
        """
        # Default to CPU count if num_workers is None
        self.num_workers = num_workers or os.cpu_count() or 4
        self.track_progress = track_progress
        self.progress_tracker = None
        
    def _create_pool(self):
        """Create the worker pool. To be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement this method")
        
    def _close_pool(self):
        """Close the worker pool. To be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement this method")
    
    def map(self, func: Callable, tasks: List[Any], chunksize: Optional[int] = None) -> List[Any]:
        """
        Apply function to each task in parallel.
        
        Args:
            func: Function to apply to each task
            tasks: List of tasks
            chunksize: Size of task chunks (optional)
            
        Returns:
            List of results
        """
        raise NotImplementedError("Subclasses must implement this method")
    
    def __enter__(self):
        """Context manager entry point."""
        self._create_pool()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit point with resource cleanup."""
        self._close_pool()


class ProcessWorkerPool(BaseWorkerPool):
    """Worker pool using multiprocessing."""
    
    def _create_pool(self):
        """Create process pool."""
        self.pool = ProcessPool(processes=self.num_workers)
        
    def _close_pool(self):
        """Close process pool."""
        if hasattr(self, 'pool'):
            self.pool.close()
            self.pool.join()
    
    def _worker_wrapper(self, args):
        """
        Wrapper for worker function to track progress.
        
        Args:
            args: Tuple of (function, task, task_id)
            
        Returns:
            Function result
        """
        func, task, task_id = args
        try:
            result = func(task)
            
            if self.progress_tracker:
                self.progress_tracker.increment()
                
            return result
        except Exception as e:
            # Catch and return exceptions to avoid process worker crashes
            if self.progress_tracker:
                self.progress_tracker.increment()
            return e
    
    def map(self, func: Callable, tasks: List[Any], chunksize: Optional[int] = None) -> List[Any]:
        """
        Apply function to each task in parallel using processes.
        
        Args:
            func: Function to apply to each task
            tasks: List of tasks
            chunksize: Size of task chunks (optional)
            
        Returns:
            List of results
        """
        if not tasks:
            return []
            
        # Create pool if not already created
        if not hasattr(self, 'pool'):
            self._create_pool()
            
        # For handling lambdas and class methods that might not pickle properly
        # We'll use a simple approach that works for most common cases
        def safe_apply(task):
            try:
                return func(task)
            except Exception as e:
                return e
        
        # Setup progress tracking
        if self.track_progress:
            self.progress_tracker = ProgressTracker(len(tasks))
            
            # Create a fixed local progress tracker to avoid sharing across processes
            total_tasks = len(tasks)
            
            # Define a self-contained function that doesn't rely on shared state
            def tracked_worker(task):
                try:
                    result = func(task)
                    # Don't use the shared progress tracker here - we'll update it in the main process
                    return result
                except Exception as e:
                    return e
            
            # Execute tasks in chunks to avoid large serialization
            all_results = []
            for i in range(0, len(tasks), 100):  # Process in smaller batches
                batch = tasks[i:i+100]
                batch_results = self.pool.map(
                    tracked_worker,
                    batch,
                    chunksize or max(1, len(batch) // (self.num_workers * 4))
                )
                all_results.extend(batch_results)
                # Update progress here in the main process
                self.progress_tracker.increment(len(batch))
                
            results = all_results
        else:
            # Execute directly without wrapper but still in batches for large task lists
            all_results = []
            for i in range(0, len(tasks), 100):  # Process in smaller batches
                batch = tasks[i:i+100]
                batch_results = self.pool.map(
                    safe_apply,
                    batch,
                    chunksize or max(1, len(batch) // (self.num_workers * 4))
                )
                all_results.extend(batch_results)
            
            results = all_results
            
        # Re-raise any exceptions that were returned
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                raise result
                
        return results


class ThreadWorkerPool(BaseWorkerPool):
    """Worker pool using threading."""
    
    def _create_pool(self):
        """Create thread pool."""
        self.pool = ThreadPoolExecutor(max_workers=self.num_workers)
        
    def _close_pool(self):
        """Close thread pool."""
        if hasattr(self, 'pool'):
            self.pool.shutdown()
    
    def map(self, func: Callable, tasks: List[Any], chunksize: Optional[int] = None) -> List[Any]:
        """
        Apply function to each task in parallel using threads.
        
        Args:
            func: Function to apply to each task
            tasks: List of tasks
            chunksize: Size of task chunks (ignored in ThreadWorkerPool)
            
        Returns:
            List of results
        """
        if not tasks:
            return []
            
        # Create pool if not already created
        if not hasattr(self, 'pool'):
            self._create_pool()
            
        # Define a safe function that catches exceptions
        def safe_apply(task):
            try:
                return func(task)
            except Exception as e:
                return e
        
        # Setup progress tracking
        if self.track_progress:
            self.progress_tracker = ProgressTracker(len(tasks))
            
            # Define wrapper function for progress tracking
            def _tracked_func(task):
                try:
                    result = func(task)
                    self.progress_tracker.increment()
                    return result
                except Exception as e:
                    self.progress_tracker.increment()
                    return e
                
            # Execute with wrapper function in smaller batches for better progress reporting
            all_results = []
            for i in range(0, len(tasks), 100):  # Process in batches
                batch = tasks[i:i+100]
                batch_results = list(self.pool.map(_tracked_func, batch))
                all_results.extend(batch_results)
            
            results = all_results
        else:
            # Execute directly without wrapper but still in batches for large task lists
            all_results = []
            for i in range(0, len(tasks), 100):
                batch = tasks[i:i+100]
                batch_results = list(self.pool.map(safe_apply, batch))
                all_results.extend(batch_results)
            
            results = all_results
            
        # Re-raise any exceptions that were returned
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                raise result
                
        return results


def create_worker_pool(pool_type: str = 'process', num_workers: Optional[int] = None, 
                      track_progress: bool = True) -> BaseWorkerPool:
    """
    Factory function to create appropriate worker pool.
    
    Args:
        pool_type: Type of pool ('process' or 'thread')
        num_workers: Number of worker processes/threads
        track_progress: Whether to track and report progress
        
    Returns:
        Worker pool instance
        
    Raises:
        ValueError: If pool_type is not recognized
    """
    if pool_type.lower() == 'process':
        return ProcessWorkerPool(num_workers, track_progress)
    elif pool_type.lower() == 'thread':
        return ThreadWorkerPool(num_workers, track_progress)
    else:
        raise ValueError(f"Unknown pool type: {pool_type}. Use 'process' or 'thread'.")


def execute_parallel(func: Callable, tasks: List[Any], 
                    num_workers: Optional[int] = None,
                    pool_type: str = 'process',
                    chunksize: Optional[int] = None,
                    track_progress: bool = True) -> List[Any]:
    """
    Convenience function to execute tasks in parallel.
    
    Args:
        func: Function to apply to each task
        tasks: List of tasks
        num_workers: Number of worker processes/threads
        pool_type: Type of pool ('process' or 'thread')
        chunksize: Size of task chunks
        track_progress: Whether to track and report progress
        
    Returns:
        List of results
    """
    with create_worker_pool(pool_type, num_workers, track_progress) as pool:
        return pool.map(func, tasks, chunksize)
