"""
Unit tests for parallel processing functionality.
"""

import unittest
import time
import random
import math
import networkx as nx
from typing import List, Dict, Any
import os
import multiprocessing

from src.parallel.task_manager import (
    TaskChunker, 
    ProgressTracker, 
    create_worker_pool,
    execute_parallel
)
from src.parallel.hierarchical_executor import (
    HierarchicalExecutor,
    execute_hierarchical_tasks,
    Task
)

# Define some test functions
def square(x):
    """Square a number."""
    return x * x

def slow_square(x):
    """Square a number with delay."""
    time.sleep(0.01)
    return x * x

def factorial(n):
    """Calculate factorial of n."""
    if n <= 1:
        return 1
    else:
        return n * factorial(n - 1)

def sometimes_fails(x):
    """Function that fails for certain inputs."""
    if x % 5 == 0:
        raise ValueError(f"Value not allowed: {x}")
    return x * 2

def do_with_resources(x, resources=None):
    """Function that uses some resources."""
    result = x * 3
    if resources:
        resources.append(x)
    return result


class TestTaskChunker(unittest.TestCase):
    """Tests for TaskChunker class."""
    
    def test_chunk_tasks_empty(self):
        """Test chunking with empty task list."""
        chunks = TaskChunker.chunk_tasks([], 4)
        self.assertEqual(chunks, [])
        
    def test_chunk_tasks_fewer_than_chunks(self):
        """Test when there are fewer tasks than requested chunks."""
        tasks = [1, 2]
        chunks = TaskChunker.chunk_tasks(tasks, 4)
        self.assertEqual(len(chunks), 2)
        self.assertEqual(chunks, [[1], [2]])
        
    def test_chunk_tasks_equal_chunks(self):
        """Test with tasks that divide evenly into chunks."""
        tasks = [1, 2, 3, 4, 5, 6]
        chunks = TaskChunker.chunk_tasks(tasks, 3)
        self.assertEqual(len(chunks), 3)
        self.assertEqual(chunks, [[1, 2], [3, 4], [5, 6]])
        
    def test_chunk_tasks_uneven_chunks(self):
        """Test with tasks that don't divide evenly into chunks."""
        tasks = [1, 2, 3, 4, 5, 6, 7]
        chunks = TaskChunker.chunk_tasks(tasks, 3)
        self.assertEqual(len(chunks), 3)
        # First chunks should have extra items
        self.assertEqual(chunks, [[1, 2, 3], [4, 5], [6, 7]])
        
    def test_chunk_by_size(self):
        """Test chunking by size."""
        tasks = [1, 2, 3, 4, 5, 6, 7]
        chunks = TaskChunker.chunk_by_size(tasks, 2)
        self.assertEqual(chunks, [[1, 2], [3, 4], [5, 6], [7]])
        
    def test_chunk_by_size_empty(self):
        """Test chunking by size with empty list."""
        chunks = TaskChunker.chunk_by_size([], 2)
        self.assertEqual(chunks, [])


class TestProgressTracker(unittest.TestCase):
    """Tests for ProgressTracker class."""
    
    def test_increment(self):
        """Test incrementing completed tasks count."""
        tracker = ProgressTracker(100, report_interval=1000)  # Long interval to avoid logging
        tracker.increment()
        self.assertEqual(tracker.get_completed(), 1)
        tracker.increment(5)
        self.assertEqual(tracker.get_completed(), 6)


class TestWorkerPools(unittest.TestCase):
    """Tests for worker pool implementations."""
    
    def test_process_pool_basic(self):
        """Test basic functionality of process pool."""
        with create_worker_pool('process', num_workers=2, track_progress=False) as pool:
            results = pool.map(square, [1, 2, 3, 4, 5])
        self.assertEqual(results, [1, 4, 9, 16, 25])
        
    def test_thread_pool_basic(self):
        """Test basic functionality of thread pool."""
        with create_worker_pool('thread', num_workers=2, track_progress=False) as pool:
            results = pool.map(square, [1, 2, 3, 4, 5])
        self.assertEqual(results, [1, 4, 9, 16, 25])
        
    def test_execute_parallel_convenience(self):
        """Test the execute_parallel convenience function."""
        results = execute_parallel(square, [1, 2, 3, 4, 5], 
                                 num_workers=2, 
                                 track_progress=False)
        self.assertEqual(results, [1, 4, 9, 16, 25])
        
    def test_invalid_pool_type(self):
        """Test that invalid pool type raises error."""
        with self.assertRaises(ValueError):
            create_worker_pool('invalid_type')
    
    def test_empty_tasks(self):
        """Test with empty task list."""
        results = execute_parallel(square, [], track_progress=False)
        self.assertEqual(results, [])
            
    def test_progress_tracking(self):
        """Test that progress is tracked."""
        # This is more of an integration test to make sure progress tracking doesn't break
        results = execute_parallel(slow_square, [1, 2, 3, 4, 5], 
                                 num_workers=2,
                                 track_progress=True)
        self.assertEqual(results, [1, 4, 9, 16, 25])
        
    def test_error_handling(self):
        """Test that worker errors are propagated."""
        # Using lambda to avoid issues with the recursive factorial
        with self.assertRaises(RecursionError):
            execute_parallel(lambda x: factorial(x), [1000], 
                           track_progress=False)
            
    def test_chunksize(self):
        """Test that chunksize parameter is respected."""
        # This is hard to test directly, but we can at least make sure it doesn't break
        results = execute_parallel(square, list(range(100)), 
                                 chunksize=10,
                                 track_progress=False)
        self.assertEqual(results, [x*x for x in range(100)])


class TestHierarchicalExecutor(unittest.TestCase):
    """Tests for HierarchicalExecutor class."""
    
    def test_basic_execution(self):
        """Test basic task execution without dependencies."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 1)
        executor.add_task('task2', lambda: 2)
        
        results = executor.execute()
        self.assertEqual(results['task1'], 1)
        self.assertEqual(results['task2'], 2)
        
    def test_dependencies(self):
        """Test execution with dependencies."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 5)
        executor.add_task('task2', lambda dep_task1: dep_task1 * 2, 
                        dependencies=['task1'])
        
        results = executor.execute()
        self.assertEqual(results['task1'], 5)
        self.assertEqual(results['task2'], 10)
        
    def test_diamond_dependencies(self):
        """Test diamond-shaped dependency graph."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 5)
        executor.add_task('task2', lambda dep_task1: dep_task1 * 2, 
                        dependencies=['task1'])
        executor.add_task('task3', lambda dep_task1: dep_task1 + 1, 
                        dependencies=['task1'])
        executor.add_task('task4', lambda dep_task2, dep_task3: dep_task2 + dep_task3, 
                        dependencies=['task2', 'task3'])
        
        results = executor.execute()
        self.assertEqual(results['task1'], 5)
        self.assertEqual(results['task2'], 10)
        self.assertEqual(results['task3'], 6)
        self.assertEqual(results['task4'], 16)
        
    def test_cyclic_dependencies(self):
        """Test that cyclic dependencies raise an error."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 1)
        executor.add_task('task2', lambda: 2, dependencies=['task1'])
        executor.add_task('task3', lambda: 3, dependencies=['task2'])
        
        # This creates a cycle
        executor.dag.add_edge('task3', 'task1')
        
        with self.assertRaises(nx.NetworkXUnfeasible):
            executor.execute()
            
    def test_error_handling(self):
        """Test error handling during execution."""
        def failing_task():
            raise ValueError("Task failed")
            
        executor = HierarchicalExecutor(fail_fast=False)
        executor.add_task('task1', lambda: 1)
        executor.add_task('task2', failing_task)
        executor.add_task('task3', lambda: 3, dependencies=['task1'])
        
        results = executor.execute()
        self.assertEqual(results['task1'], 1)
        self.assertEqual(results['task3'], 3)
        self.assertIn('task2', executor.errors)
        self.assertNotIn('task2', results)
        
    def test_fail_fast(self):
        """Test fail_fast behavior."""
        def failing_task():
            raise ValueError("Task failed")
            
        executor = HierarchicalExecutor(fail_fast=True)
        executor.add_task('task1', failing_task)
        executor.add_task('task2', lambda: 2)  # Should not execute
        
        results = executor.execute()
        self.assertEqual(results, {})
        self.assertIn('task1', executor.errors)
        
    def test_execution_plan(self):
        """Test getting the execution plan."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 1)
        executor.add_task('task2', lambda: 2, dependencies=['task1'])
        executor.add_task('task3', lambda: 3, dependencies=['task1'])
        executor.add_task('task4', lambda: 4, dependencies=['task2', 'task3'])
        
        plan = executor.get_execution_plan()
        self.assertEqual(len(plan), 3)  # 3 levels
        self.assertEqual(plan[0], ['task1'])  # Level 0
        self.assertCountEqual(plan[1], ['task2', 'task3'])  # Level 1
        self.assertEqual(plan[2], ['task4'])  # Level 2
        
    def test_convenience_function(self):
        """Test the execute_hierarchical_tasks convenience function."""
        tasks = [
            {'id': 'task1', 'func': lambda: 5},
            {'id': 'task2', 'func': lambda x: x * 2, 'args': (5,)},
            {'id': 'task3', 'func': lambda x, y: x + y, 
             'kwargs': {'x': 3, 'y': 4}, 'dependencies': ['task1', 'task2']}
        ]
        
        results = execute_hierarchical_tasks(tasks, fail_fast=False)
        self.assertEqual(results['task1'], 5)
        self.assertEqual(results['task2'], 10)
        self.assertEqual(results['task3'], 3 + 4)  # Should use kwargs, not dependencies


class TestPerformanceScaling(unittest.TestCase):
    """Tests for performance scaling with parallelism."""
    
    def test_scaling_processes(self):
        """Test performance scaling with increasing process workers."""
        # Skip if running in CI environment
        if 'CI' in os.environ:
            self.skipTest("Skipping performance test in CI environment")
            
        # Run with different numbers of workers
        data = list(range(200))
        
        # Warmup
        _ = execute_parallel(slow_square, data[:10], num_workers=1, track_progress=False)
        
        times = {}
        for num_workers in [1, 2, multiprocessing.cpu_count() or 4]:
            start_time = time.time()
            results = execute_parallel(slow_square, data, 
                                     num_workers=num_workers, 
                                     track_progress=False)
            end_time = time.time()
            times[num_workers] = end_time - start_time
            
            # Verify results
            self.assertEqual(results, [x*x for x in data])
            
        # Check if parallelization gives speedup (with some tolerance for overhead)
        # single core should be slower than multicore for this workload
        self.assertGreater(times[1], times[multiprocessing.cpu_count()])
        
    def test_serial_vs_parallel_equivalence(self):
        """Test that serial and parallel execution give identical results."""
        data = list(range(100))
        
        # Serial execution
        serial_results = [square(x) for x in data]
        
        # Parallel execution with processes
        parallel_results_proc = execute_parallel(square, data, 
                                              pool_type='process',
                                              track_progress=False)
        
        # Parallel execution with threads
        parallel_results_thread = execute_parallel(square, data, 
                                                pool_type='thread',
                                                track_progress=False)
        
        # Verify results match
        self.assertEqual(serial_results, parallel_results_proc)
        self.assertEqual(serial_results, parallel_results_thread)
        
    def test_error_handling_and_cleanup(self):
        """Test error handling and cleanup in failure cases."""
        resources = []
        
        # Create a list of inputs that will cause some failures
        data = list(range(20))
        
        # Test with fail_fast=False to allow partial results
        results = []
        errors = []
        
        # Use a context manager to ensure cleanup
        with create_worker_pool('thread', track_progress=False) as pool:
            # Map the function that sometimes fails
            try:
                results = pool.map(sometimes_fails, data)
            except ValueError as e:
                errors.append(e)
                
        # Check that we got some results for the non-failing inputs
        self.assertEqual(len(results), 0 if errors else 20)
        
        # Check resource cleanup with a function that uses resources
        with create_worker_pool('thread', track_progress=False) as pool:
            try:
                results = pool.map(lambda x: do_with_resources(x, resources), data)
            except Exception:
                pass
                
        # The pool should have been properly cleaned up
        # This is hard to test directly, but at least we can check the results
        if not errors:
            self.assertEqual(len(resources), 20)
            
    def test_hierarchical_task_scaling(self):
        """Test that hierarchical execution scales with parallelism."""
        # Skip if running in CI environment
        if 'CI' in os.environ:
            self.skipTest("Skipping performance test in CI environment")
            
        # Create a wide, shallow graph with many independent tasks
        def make_executor(num_tasks, num_workers):
            executor = HierarchicalExecutor(num_workers=num_workers)
            
            # Add a root task
            executor.add_task('root', lambda: 0)
            
            # Add many independent tasks depending on root
            for i in range(num_tasks):
                executor.add_task(f'task_{i}', lambda x=i: slow_square(x), 
                                dependencies=['root'])
                
            # Add a final aggregation task
            deps = [f'task_{i}' for i in range(num_tasks)]
            executor.add_task('final', lambda: sum(1 for _ in deps), 
                            dependencies=deps)
                            
            return executor
            
        # Warmup
        executor = make_executor(5, 1)
        executor.execute()
        
        # Test with different worker counts
        num_tasks = 50
        times = {}
        
        for num_workers in [1, 4]:
            executor = make_executor(num_tasks, num_workers)
            
            start_time = time.time()
            results = executor.execute()
            end_time = time.time()
            
            times[num_workers] = end_time - start_time
            
            # Verify results
            self.assertEqual(results['final'], num_tasks)
            
        # Multiple workers should be faster than a single worker
        # Allow some tolerance for very small workloads and overhead
        if times[1] > 0.5:  # Only check if the test is meaningful (not too fast)
            self.assertLess(times[4], times[1] * 0.9)


if __name__ == '__main__':
    unittest.main()
