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
import logging # Import logging to use assertLogs
import sys # Import sys for stream handler

# Ensure src is in path (relative path from tests/ to src/)
script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(script_dir, '..', 'src'))
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

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
        # Limit recursion depth for testing purposes
        if n > 500:
             raise RecursionError("Maximum recursion depth exceeded (simulated)")
        return n * factorial(n - 1)

# Top-level function for pickling compatibility with multiprocessing
def _test_factorial(n):
    """Helper function to call factorial for testing."""
    return factorial(n)

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

    def setUp(self):
        """Set up basic logging for pool tests."""
        # Store original handlers and level to restore later
        self.original_handlers = logging.getLogger().handlers[:]
        self.original_level = logging.getLogger().level
        # Configure basic logging to capture messages during tests
        logging.basicConfig(level=logging.DEBUG, stream=sys.stderr, force=True)

    def tearDown(self):
        """Clean up logging configuration."""
        root_logger = logging.getLogger()
        for handler in root_logger.handlers[:]:
             root_logger.removeHandler(handler)
             handler.close()
        for handler in self.original_handlers:
             root_logger.addHandler(handler)
        root_logger.setLevel(self.original_level)

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
        # Store original inputs to verify output order
        inputs = [1, 2, 3, 4, 5]
        results = execute_parallel(square, inputs,
                                 num_workers=2,
                                 track_progress=False)
        # Create expected outputs in the same order as inputs
        expected = [square(x) for x in inputs]
        self.assertEqual(results, expected)

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
        # Use root logger since tests may be run in various environments
        # where logger names might differ
        with self.assertLogs(level='INFO') as cm:
            results = execute_parallel(slow_square, [1, 2, 3, 4, 5],
                                     num_workers=2,
                                     track_progress=True)
        expected = [slow_square(x) for x in [1, 2, 3, 4, 5]]
        self.assertEqual(results, expected)
        # Check that at least one progress message was logged
        self.assertTrue(any("Progress:" in log for log in cm.output), f"No progress logs found in {cm.output}")


    def test_error_handling(self):
        """Test that worker errors are propagated."""
        # Use top-level function for process pool compatibility
        with self.assertRaises(RecursionError):
            # Use process pool to ensure error pickling works
            execute_parallel(_test_factorial, [1000], pool_type='process',
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

    def setUp(self):
        """Set up basic logging for HierarchicalExecutor tests."""
        # Store original handlers and level to restore later
        self.original_handlers = logging.getLogger().handlers[:]
        self.original_level = logging.getLogger().level
        # Configure basic logging to capture messages during tests
        logging.basicConfig(level=logging.DEBUG, stream=sys.stderr, force=True)

    def tearDown(self):
        """Clean up logging configuration."""
        root_logger = logging.getLogger()
        for handler in root_logger.handlers[:]:
             root_logger.removeHandler(handler)
             handler.close()
        for handler in self.original_handlers:
             root_logger.addHandler(handler)
        root_logger.setLevel(self.original_level)

    def test_basic_execution(self):
        """Test basic task execution without dependencies."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 1)
        executor.add_task('task2', lambda: 2)

        results = executor.execute()
        self.assertEqual(results.get('task1'), 1) # Use .get() for safer access
        self.assertEqual(results.get('task2'), 2)
        self.assertEqual(len(results), 2) # Ensure only expected results are present

    def test_dependencies(self):
        """Test execution with dependencies."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 5)
        # Lambda now expects dependency result positionally *after* own args (if any)
        # Since task2 has no args, dep_task1_result is the first arg.
        executor.add_task('task2', lambda dep_task1_result: dep_task1_result * 2,
                        dependencies=['task1'])

        results = executor.execute()
        self.assertEqual(results.get('task1'), 5)
        self.assertEqual(results.get('task2'), 10)
        self.assertEqual(len(results), 2)

    def test_diamond_dependencies(self):
        """Test diamond-shaped dependency graph."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 5)
        # Lambdas expect dependency results positionally
        # task2 depends on task1
        executor.add_task('task2', lambda dep_task1_result: dep_task1_result * 2,
                        dependencies=['task1'])
        # task3 depends on task1
        executor.add_task('task3', lambda dep_task1_result: dep_task1_result + 1,
                        dependencies=['task1'])
        # task4 depends on task2 and task3 (results passed positionally in order of dependencies list)
        executor.add_task('task4', lambda dep_task2_result, dep_task3_result: dep_task2_result + dep_task3_result,
                        dependencies=['task2', 'task3'])

        results = executor.execute()
        self.assertEqual(results.get('task1'), 5)
        self.assertEqual(results.get('task2'), 10)
        self.assertEqual(results.get('task3'), 6)
        self.assertEqual(results.get('task4'), 16)
        self.assertEqual(len(results), 4)

    def test_cyclic_dependencies(self):
        """Test that cyclic dependencies raise an error during execution."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 1, dependencies=['task3']) # Create cycle: 1 -> 2 -> 3 -> 1
        executor.add_task('task2', lambda dep1: 2, dependencies=['task1'])
        executor.add_task('task3', lambda dep2: 3, dependencies=['task2'])

        # The cycle check happens in execute() or get_execution_plan()
        with self.assertRaises(nx.NetworkXUnfeasible):
            executor.execute()
        with self.assertRaises(nx.NetworkXUnfeasible):
            executor.get_execution_plan()


    def test_error_handling(self):
        """Test error handling during execution without fail_fast."""
        def failing_task():
            raise ValueError("Task failed intentionally")

        executor = HierarchicalExecutor(fail_fast=False)
        executor.add_task('task1', lambda: 1)
        executor.add_task('task2', failing_task) # This task will fail
        # Task 3 depends on task 1 (should run)
        executor.add_task('task3', lambda dep1_res: dep1_res + 2, dependencies=['task1'])
        # Task 4 depends on task 2 (should fail due to dependency)
        executor.add_task('task4', lambda dep2_res: dep2_res + 5, dependencies=['task2'])
        # Task 5 depends on task 3 (should run)
        executor.add_task('task5', lambda dep3_res: dep3_res * 3, dependencies=['task3'])

        # Capture expected ERROR log for task2 and WARNING logs for task4
        # Target the specific logger used by the executor
        with self.assertLogs('src.parallel.hierarchical_executor', level='WARNING') as cm: # Capture WARNING and above
            results = executor.execute()

        # Check successful tasks
        self.assertEqual(results.get('task1'), 1)
        self.assertEqual(results.get('task3'), 3) # 1 + 2 = 3
        self.assertEqual(results.get('task5'), 9) # 3 * 3 = 9

        # Check failed tasks
        self.assertIn('task2', executor.errors)
        self.assertIsInstance(executor.errors['task2'], ValueError)
        self.assertNotIn('task2', results)

        self.assertIn('task4', executor.errors) # Should fail due to dependency
        self.assertIsInstance(executor.errors['task4'], RuntimeError) # Dependency failure error
        self.assertNotIn('task4', results)

        # Verify the expected log messages were captured
        log_output = "\n".join(cm.output)
        self.assertIn("ERROR:src.parallel.hierarchical_executor:Task execution failed: task2", log_output)
        self.assertIn("ValueError: Task failed intentionally", log_output)
        self.assertIn("WARNING:src.parallel.hierarchical_executor:Task failed: task2", log_output)
        self.assertIn("WARNING:src.parallel.hierarchical_executor:Marking task task4 as failed due to upstream dependency failure on task2", log_output)


    def test_fail_fast(self):
        """Test fail_fast behavior, capturing expected error logs."""
        def failing_task():
            time.sleep(0.01) # Ensure it doesn't finish instantly
            raise ValueError("Task failed intentionally for fail_fast")

        def slow_task():
            time.sleep(0.1)
            return 2

        executor = HierarchicalExecutor(fail_fast=True, num_workers=2)
        executor.add_task('task1', failing_task) # Will fail
        executor.add_task('task2', slow_task)  # Should be cancelled or not complete
        executor.add_task('task3', lambda: 3) # Might start but should be cancelled
        executor.add_task('task4', lambda dep3_res: dep3_res + 1, dependencies=['task3']) # Should not run

        # Use assertLogs to capture the expected ERROR/WARNING messages
        with self.assertLogs('src.parallel.hierarchical_executor', level='WARNING') as cm:
            results = executor.execute()

        # Check that the failing task's error is recorded
        self.assertIn('task1', executor.errors)
        self.assertIsInstance(executor.errors['task1'], ValueError)
        self.assertNotIn('task1', results)

        # Check that other tasks were either not run or marked as failed/cancelled
        self.assertNotIn('task2', results) # slow_task likely didn't finish or was cancelled
        self.assertNotIn('task3', results)
        self.assertNotIn('task4', results)

        # Check errors dictionary for cancelled tasks (they should have RuntimeError)
        self.assertIn('task2', executor.errors) # May or may not be present depending on timing
        self.assertIn('task3', executor.errors)
        self.assertIn('task4', executor.errors)
        if 'task2' in executor.errors: self.assertIsInstance(executor.errors['task2'], (RuntimeError, ValueError)) # Could be cancelled or fail itself if started
        if 'task3' in executor.errors: self.assertIsInstance(executor.errors['task3'], RuntimeError)
        if 'task4' in executor.errors: self.assertIsInstance(executor.errors['task4'], RuntimeError)


        # Verify that the expected error messages were logged
        log_output = "\n".join(cm.output)
        # Error from the task itself
        self.assertIn("ERROR:src.parallel.hierarchical_executor:Task execution failed: task1", log_output)
        self.assertIn("ValueError: Task failed intentionally for fail_fast", log_output)
        # Warning about the failure
        self.assertIn("WARNING:src.parallel.hierarchical_executor:Task failed: task1", log_output)
        # Error about fail_fast stopping execution
        self.assertIn("ERROR:src.parallel.hierarchical_executor:Fail fast enabled. Stopping execution due to task failure: task1", log_output)
        # Warnings about marking other tasks as cancelled/failed
        self.assertTrue(any("Marking task task" in log and "as cancelled due to fail_fast" in log for log in cm.output),
                        f"Expected cancellation logs not found in {log_output}")


    def test_execution_plan(self):
        """Test getting the execution plan."""
        executor = HierarchicalExecutor()
        executor.add_task('task1', lambda: 1)
        executor.add_task('task2', lambda dep1: 2, dependencies=['task1'])
        executor.add_task('task3', lambda dep1: 3, dependencies=['task1'])
        executor.add_task('task4', lambda dep2, dep3: 4, dependencies=['task2', 'task3'])

        plan = executor.get_execution_plan()
        self.assertEqual(len(plan), 3)  # 3 levels
        self.assertCountEqual(plan[0], ['task1'])  # Level 0 (use assertCountEqual for order independence)
        self.assertCountEqual(plan[1], ['task2', 'task3'])  # Level 1
        self.assertCountEqual(plan[2], ['task4'])  # Level 2

    def test_convenience_function(self):
        """Test the execute_hierarchical_tasks convenience function."""
        tasks = [
            {'id': 'task1', 'func': lambda: 5},
            # Task2 uses its own args (x=10), no dependencies
            {'id': 'task2', 'func': lambda x: x * 2, 'args': (10,)},
            # Task3 depends on task1 and task2, uses their results positionally,
            # and also uses its own kwargs.
            # Lambda signature: (dep1_res, dep2_res, z=kwarg)
            {'id': 'task3',
             'func': lambda dep1_res, dep2_res, z: dep1_res + dep2_res + z,
             'kwargs': {'z': 1},
             'dependencies': ['task1', 'task2']}
        ]

        results = execute_hierarchical_tasks(tasks, fail_fast=False)

        # Check results: task1=5, task2=20 (from 10*2), task3=26 (from 5 + 20 + 1)
        self.assertEqual(results.get('task1'), 5)
        self.assertEqual(results.get('task2'), 20)
        self.assertEqual(results.get('task3'), 26)
        self.assertEqual(len(results), 3)


class TestPerformanceScaling(unittest.TestCase):
    """Tests for performance scaling with parallelism."""

    def setUp(self):
        """Set up basic logging for performance tests."""
        self.original_handlers = logging.getLogger().handlers[:]
        self.original_level = logging.getLogger().level
        logging.basicConfig(level=logging.INFO, stream=sys.stderr, force=True) # Use INFO for performance tests

    def tearDown(self):
        """Clean up logging configuration."""
        root_logger = logging.getLogger()
        for handler in root_logger.handlers[:]:
             root_logger.removeHandler(handler)
             handler.close()
        for handler in self.original_handlers:
             root_logger.addHandler(handler)
        root_logger.setLevel(self.original_level)

    @unittest.skipIf('CI' in os.environ, "Skipping performance test in CI environment")
    def test_scaling_processes(self):
        """Test performance scaling with increasing process workers."""

        # Run with different numbers of workers
        data = list(range(200))

        # Warmup
        _ = execute_parallel(slow_square, data[:10], num_workers=1, track_progress=False)

        times = {}
        cpu_count = multiprocessing.cpu_count() or 4
        workers_to_test = sorted(list(set([1, 2, cpu_count]))) # Ensure unique, sorted workers [1, 2, cpu_count]

        for num_workers in workers_to_test:
            start_time = time.time()
            results = execute_parallel(slow_square, data,
                                     num_workers=num_workers,
                                     pool_type='process', # Explicitly use process
                                     track_progress=False)
            end_time = time.time()
            times[num_workers] = end_time - start_time

            # Verify results
            self.assertEqual(results, [x*x for x in data])

        # Check if parallelization gives speedup (with some tolerance for overhead)
        # single core should be slower than multicore for this workload
        print(f"\nScaling times (processes): {times}") # Optional: print times for debugging
        if 1 in times and cpu_count in times and times[1] > 0.1: # Only assert if runtime is significant enough
             # Expect multicore to be faster than single core
             self.assertLess(times[cpu_count], times[1], f"Expected speedup: {cpu_count} workers ({times[cpu_count]:.3f}s) should be faster than 1 worker ({times[1]:.3f}s)")
             # Expect 2 workers to be faster than 1 worker (if cpu_count > 1)
             if 2 in times and cpu_count > 1:
                  self.assertLess(times[2], times[1], f"Expected speedup: 2 workers ({times[2]:.3f}s) should be faster than 1 worker ({times[1]:.3f}s)")


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
        """Test error handling and cleanup in failure cases for standard pools."""
        # Use a mutable object like a list for resources if needed across threads/processes
        # For simple tests, direct check might suffice

        # Create a list of inputs that will cause some failures
        data = list(range(20)) # Will fail for 0, 5, 10, 15

        # Test with fail_fast=False (default for execute_parallel)
        caught_exception = None

        # Use a context manager to ensure cleanup
        try:
            # Use 'thread' pool to avoid pickling issues with sometimes_fails if run inside method
            # Capture expected ERROR logs from the pool map without printing them
            # Note: execute_parallel itself doesn't log errors from the mapped function,
            # it relies on the pool implementation to raise them.
            # So, we don't use assertLogs here, but expect an exception.
            _ = execute_parallel(sometimes_fails, data, pool_type='thread', track_progress=False)

        except ValueError as e:
             caught_exception = e
        except Exception as e:
             self.fail(f"Expected ValueError, but got {type(e).__name__}: {e}")


        # Check that an exception was indeed raised (the first one encountered)
        self.assertIsInstance(caught_exception, ValueError)
        # The first failure should be for x=0
        self.assertIn("Value not allowed: 0", str(caught_exception))


        # Test resource cleanup (less direct, more about ensuring pool closes)
        resources_list = [] # Use a list accessible in this scope
        try:
             # Use lambda to pass the list correctly
             _ = execute_parallel(lambda x: do_with_resources(x, resources=resources_list),
                                  data, pool_type='thread', track_progress=False)
        except Exception:
             pass # Ignore exceptions for this part, focus on cleanup

        # Check if resources were potentially populated (best effort check)
        # This might be empty if an error occurred early
        # A better check might be to ensure the pool context manager exited cleanly
        # If no error occurred, all items should be processed
        if caught_exception is None: # If the previous test didn't catch an error (e.g., sometimes_fails changed)
            self.assertEqual(len(resources_list), len(data))
        else:
            # If an error occurred, we can't guarantee how many tasks completed
            self.assertTrue(len(resources_list) <= len(data))


    @unittest.skipIf('CI' in os.environ, "Skipping performance test in CI environment")
    def test_hierarchical_task_scaling(self):
        """Test that hierarchical execution scales with parallelism."""

        # Create a wide, shallow graph with many independent tasks
        def make_executor(num_tasks, num_workers):
            executor = HierarchicalExecutor(num_workers=num_workers)

            # Add a root task
            executor.add_task('root', lambda: 0)

            # Add many independent tasks depending on root
            for i in range(num_tasks):
                # Lambda now expects root result positionally (though it's unused here)
                # Use default argument for i to capture its value at definition time
                executor.add_task(f'task_{i}', lambda root_res, x=i: slow_square(x),
                                dependencies=['root'])

            # Add a final aggregation task
            deps = [f'task_{i}' for i in range(num_tasks)]
            # Lambda now expects positional args for all dependencies
            # The number of args must match len(deps)
            executor.add_task('final', lambda *dep_results: sum(dep_results), # Example: sum results
                            dependencies=deps)

            return executor

        # Warmup
        executor_warmup = make_executor(5, 1)
        executor_warmup.execute()

        # Test with different worker counts
        num_tasks = 50
        times = {}
        workers_to_test = sorted(list(set([1, 4, multiprocessing.cpu_count() or 4]))) # Test 1, 4, max

        expected_final_result = sum(i*i for i in range(num_tasks))

        for num_workers in workers_to_test:
            executor = make_executor(num_tasks, num_workers)

            start_time = time.time()
            results = executor.execute()
            end_time = time.time()

            times[num_workers] = end_time - start_time

            # Verify results
            self.assertEqual(results.get('final'), expected_final_result)
            # Verify one intermediate task result as a sanity check
            self.assertEqual(results.get('task_10'), 100) # slow_square(10)

        # Multiple workers should be faster than a single worker
        # Allow some tolerance for very small workloads and overhead
        print(f"\nScaling times (hierarchical): {times}") # Optional: print times
        if 1 in times and workers_to_test[-1] in times and times[1] > 0.1:  # Only check if the test is meaningful
            self.assertLess(times[workers_to_test[-1]], times[1], f"Expected speedup: {workers_to_test[-1]} workers ({times[workers_to_test[-1]]:.3f}s) should be faster than 1 worker ({times[1]:.3f}s)")
            if 4 in times and workers_to_test[-1] >= 4:
                 self.assertLess(times[4], times[1], f"Expected speedup: 4 workers ({times[4]:.3f}s) should be faster than 1 worker ({times[1]:.3f}s)")


if __name__ == '__main__':
    # Configure logging level for tests if needed, e.g., to see INFO messages
    logging.basicConfig(level=logging.INFO, stream=sys.stderr, force=True)
    unittest.main()
