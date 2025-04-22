"""
Hierarchical executor module.

Implements parallel execution that respects hierarchical dependencies
using a directed acyclic graph to manage execution order.
"""

import logging
import time
import os  # Added import
import networkx as nx
from typing import Dict, List, Set, Callable, Any, Optional, Tuple, Union
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed, Future

from src.parallel.task_manager import create_worker_pool

# Use the logger for this specific module
logger = logging.getLogger(__name__)

@dataclass
class Task:
    """Represents a task in the hierarchical execution graph."""

    id: str
    func: Callable
    args: Tuple = ()
    kwargs: Dict = None
    dependencies: List[str] = None
    result: Any = None
    status: str = 'pending'  # 'pending', 'running', 'completed', 'failed'
    error: Optional[Exception] = None # Use Optional for clarity

    def __post_init__(self):
        """Initialize default values."""
        if self.kwargs is None:
            self.kwargs = {}
        if self.dependencies is None:
            self.dependencies = []

    def execute(self, dependency_args: Tuple = ()) -> Any:
        """
        Execute the task.

        Args:
            dependency_args: Tuple containing results from dependency tasks,
                             in the order specified by self.dependencies.

        Returns:
            Result of the task

        Raises:
            Any exception raised during task execution
        """
        self.status = 'running'
        logger.debug(f"Starting execution of task: {self.id}")
        try:
            # Combine task's own args with dependency results
            final_args = self.args + dependency_args

            # Execute the function with combined arguments
            self.result = self.func(*final_args, **self.kwargs)
            self.status = 'completed'
            logger.debug(f"Task completed successfully: {self.id}")
            return self.result

        except Exception as e:
            self.status = 'failed'
            self.error = e
            # Log the error here when it occurs during execution
            logger.error(f"Task execution failed: {self.id} - {type(e).__name__}: {e}", exc_info=True)
            raise # Re-raise the exception to be caught by the executor


class HierarchicalExecutor:
    """
    Executes tasks with respect to their hierarchical dependencies.

    Uses a directed acyclic graph (DAG) to track and manage dependencies
    between tasks, and executes independent tasks in parallel.
    """

    def __init__(self,
                num_workers: Optional[int] = None,
                fail_fast: bool = False):
        """
        Initialize hierarchical executor.

        Args:
            num_workers: Number of worker threads
            fail_fast: Whether to stop all execution on first failure
        """
        self.num_workers = num_workers
        self.fail_fast = fail_fast
        self.tasks: Dict[str, Task] = {}
        self.dag = nx.DiGraph()
        self.results: Dict[str, Any] = {}
        self.errors: Dict[str, Exception] = {}

    def add_task(self,
                task_id: str,
                func: Callable,
                args: Tuple = (),
                kwargs: Dict = None,
                dependencies: List[str] = None) -> None:
        """
        Add a task to the execution graph.

        Args:
            task_id: Unique identifier for the task
            func: Function to execute
            args: Positional arguments for the function
            kwargs: Keyword arguments for the function
            dependencies: List of task IDs this task depends on

        Raises:
            ValueError: If task_id already exists or if a dependency is unknown
        """
        if task_id in self.tasks:
            raise ValueError(f"Task ID already exists: {task_id}")

        # Create task object
        task = Task(id=task_id, func=func, args=args, kwargs=kwargs, dependencies=dependencies or [])
        self.tasks[task_id] = task

        # Add to DAG
        self.dag.add_node(task_id)

        # Add dependency edges
        for dep_id in task.dependencies:
            # Check if dependency exists *before* adding edge
            if dep_id not in self.tasks:
                 # If adding tasks out of order, allow adding node first
                 if dep_id not in self.dag:
                      self.dag.add_node(dep_id) # Add dependency node placeholder
                 # Check again after potentially adding node
                 # This logic might be complex if tasks are added totally out of order.
                 # A safer approach is to require dependencies to be added first,
                 # or validate the graph fully before execution.
                 # For now, assume dependencies are added reasonably or exist.
                 # Let's raise error if dependency task object isn't known *yet*.
                 # if dep_id not in self.tasks:
                 #     raise ValueError(f"Unknown dependency task object: {dep_id} for task {task_id}")

            self.dag.add_edge(dep_id, task_id)

    def execute(self) -> Dict[str, Any]:
        """
        Execute all tasks respecting dependencies.

        Returns:
            Dictionary mapping task IDs to results

        Raises:
            nx.NetworkXUnfeasible: If the dependency graph has cycles
            RuntimeError: If execution fails and fail_fast is True
        """
        # Validate all nodes in DAG correspond to tasks
        missing_tasks = [node for node in self.dag.nodes() if node not in self.tasks]
        if missing_tasks:
            raise ValueError(f"Graph contains nodes without corresponding task objects: {missing_tasks}")

        # Check for cycles
        try:
            if not nx.is_directed_acyclic_graph(self.dag):
                cycles = list(nx.simple_cycles(self.dag))
                raise nx.NetworkXUnfeasible(f"Dependency graph contains cycles: {cycles}")
        except nx.NetworkXNotImplemented:
             # is_directed_acyclic_graph might raise this for MultiDiGraph, but we use DiGraph
             logger.warning("Could not check for cycles in the dependency graph.")


        # Clear previous results
        self.results = {}
        self.errors = {}

        # Get topological order (ensures dependencies are executed first)
        try:
            # Use a generator for potentially large graphs
            execution_plan_generator = nx.topological_sort(self.dag)
        except nx.NetworkXUnfeasible as e:
             # This catches cycles identified by topological_sort
             raise nx.NetworkXUnfeasible(f"Dependency graph contains cycles. Cannot execute. Details: {e}") from e


        # Set up worker pool
        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            # Keep track of running futures: task_id -> Future
            futures: Dict[str, Future] = {}
            # Keep track of tasks ready to be submitted but waiting for worker
            ready_to_submit: List[str] = []
            # Keep track of tasks whose dependencies are met
            completed_dependencies: Dict[str, Set[str]] = {task_id: set() for task_id in self.tasks}
            # Tasks with no dependencies are ready immediately
            tasks_in_progress_or_done = set()

            # Find initial ready tasks (no dependencies)
            for task_id in self.tasks:
                 if not self.tasks[task_id].dependencies:
                      ready_to_submit.append(task_id)
                      tasks_in_progress_or_done.add(task_id)


            while tasks_in_progress_or_done != set(self.tasks.keys()):
                # Submit ready tasks up to the number of workers
                while ready_to_submit and len(futures) < (self.num_workers or os.cpu_count() or 4):
                    task_id = ready_to_submit.pop(0)
                    task = self.tasks[task_id]
                    logger.debug(f"Submitting task {task_id} for execution.")
                    # Dependency args should be empty for initial tasks
                    future = executor.submit(task.execute, ())
                    futures[task_id] = future

                # Wait for at least one future to complete
                if not futures:
                     # This might happen if there's a gap in the graph or an error state
                     logger.warning("No active futures, but not all tasks are done. Checking state.")
                     # Check if remaining tasks have unmet dependencies
                     all_task_ids = set(self.tasks.keys())
                     remaining_task_ids = all_task_ids - tasks_in_progress_or_done
                     if not remaining_task_ids: break # All tasks accounted for

                     can_run = False
                     for rem_task_id in remaining_task_ids:
                          task = self.tasks[rem_task_id]
                          deps_met = all(dep in self.results for dep in task.dependencies)
                          if deps_met:
                               can_run = True
                               if rem_task_id not in ready_to_submit:
                                    ready_to_submit.append(rem_task_id)
                                    tasks_in_progress_or_done.add(rem_task_id)
                     if not can_run:
                          logger.error(f"Stall detected. Remaining tasks {remaining_task_ids} have unmet dependencies or failed dependencies.")
                          # Check for failed dependencies
                          for rem_task_id in remaining_task_ids:
                               task = self.tasks[rem_task_id]
                               failed_deps = [dep for dep in task.dependencies if dep in self.errors]
                               if failed_deps:
                                    logger.error(f"Task {rem_task_id} cannot run due to failed dependencies: {failed_deps}")
                                    # Mark this task as failed due to dependency
                                    self.errors[rem_task_id] = RuntimeError(f"Dependency failure: {failed_deps}")
                                    tasks_in_progress_or_done.add(rem_task_id) # Mark as 'done' (failed)

                          # If still stalled after checking failed deps, break
                          if not ready_to_submit:
                               break
                     continue # Go back to submit loop if new tasks became ready


                # Process completed futures using as_completed for efficiency
                try:
                    # Wait for the next future to complete
                    done_futures = {f for f in as_completed(futures.values(), timeout=0.1) if f.done()}
                except TimeoutError:
                     continue # No futures completed in timeout, continue waiting

                if not done_futures: continue

                for future in done_futures:
                    # Find the task_id associated with this future
                    task_id = None
                    for tid, f in futures.items():
                        if f == future:
                            task_id = tid
                            break
                    if task_id is None: continue # Should not happen

                    del futures[task_id] # Remove from active futures

                    try:
                        result = future.result() # Get result or raise exception
                        self.results[task_id] = result
                        logger.debug(f"Task completed: {task_id}")

                        # Find successor tasks and check if they are now ready
                        for successor_id in self.dag.successors(task_id):
                             if successor_id not in tasks_in_progress_or_done:
                                  successor_task = self.tasks[successor_id]
                                  # Check if all dependencies of the successor are now met
                                  deps_met = all(dep in self.results for dep in successor_task.dependencies)
                                  if deps_met:
                                       if successor_id not in ready_to_submit:
                                            ready_to_submit.append(successor_id)
                                            tasks_in_progress_or_done.add(successor_id)

                    except Exception as e:
                        self.errors[task_id] = e
                        # Error already logged inside task.execute
                        # logger.error(f"Task failed: {task_id} - {type(e).__name__}: {e}") # Redundant log?
                        tasks_in_progress_or_done.add(task_id) # Mark as done (failed)

                        if self.fail_fast:
                            logger.error(f"Fail fast enabled. Stopping execution due to task failure: {task_id}")
                            # Cancel remaining futures
                            for f in futures.values():
                                f.cancel()
                            # Mark remaining ready tasks as cancelled/failed?
                            for rem_task_id in ready_to_submit:
                                 self.errors[rem_task_id] = RuntimeError(f"Execution stopped due to fail_fast on task {task_id}")
                                 tasks_in_progress_or_done.add(rem_task_id)
                            # Mark tasks still waiting for dependencies as failed
                            all_task_ids = set(self.tasks.keys())
                            remaining_waiting_tasks = all_task_ids - tasks_in_progress_or_done
                            for wait_task_id in remaining_waiting_tasks:
                                 self.errors[wait_task_id] = RuntimeError(f"Execution stopped due to fail_fast on task {task_id}")
                                 tasks_in_progress_or_done.add(wait_task_id)

                            return self.results # Return immediately

                        # If not fail_fast, mark dependent tasks as failed
                        queue = list(self.dag.successors(task_id))
                        visited = {task_id}
                        while queue:
                             current_dep_fail_id = queue.pop(0)
                             if current_dep_fail_id not in visited:
                                  visited.add(current_dep_fail_id)
                                  if current_dep_fail_id not in self.errors: # Don't overwrite existing errors
                                       logger.warning(f"Marking task {current_dep_fail_id} as failed due to dependency failure on {task_id}")
                                       self.errors[current_dep_fail_id] = RuntimeError(f"Dependency failure: {task_id}")
                                       tasks_in_progress_or_done.add(current_dep_fail_id)
                                  queue.extend(list(self.dag.successors(current_dep_fail_id)))


        # Final check for any tasks that might not have been marked done
        all_task_ids = set(self.tasks.keys())
        missing_from_done = all_task_ids - tasks_in_progress_or_done
        if missing_from_done:
             logger.warning(f"Some tasks were not marked as done/failed: {missing_from_done}. Marking as error.")
             for task_id in missing_from_done:
                  if task_id not in self.errors:
                       self.errors[task_id] = RuntimeError("Task did not complete execution.")


        return self.results


    def get_execution_plan(self) -> List[List[str]]:
        """
        Get the execution plan showing which tasks can execute in parallel.

        Returns:
            List of lists, where each inner list contains tasks that can
            execute in parallel at the same level.
        """
        # Check for cycles first
        try:
            if not nx.is_directed_acyclic_graph(self.dag):
                cycles = list(nx.simple_cycles(self.dag))
                raise nx.NetworkXUnfeasible(f"Dependency graph contains cycles: {cycles}")
        except nx.NetworkXNotImplemented:
             logger.warning("Could not check for cycles in the dependency graph.")
             # Proceed cautiously

        # Group nodes by generation (distance from root nodes)
        generations = list(nx.topological_generations(self.dag))
        return [list(gen) for gen in generations]


def execute_hierarchical_tasks(tasks_definitions: List[Dict],
                              num_workers: Optional[int] = None,
                              fail_fast: bool = False) -> Dict[str, Any]:
    """
    Convenience function to execute hierarchical tasks.

    Args:
        tasks_definitions: List of task definitions, each as a dictionary with:
                         - id: Task ID
                         - func: Function to execute
                         - args: Positional arguments (optional)
                         - kwargs: Keyword arguments (optional)
                         - dependencies: List of dependency task IDs (optional)
        num_workers: Number of worker threads
        fail_fast: Whether to stop all execution on first failure

    Returns:
        Dictionary mapping task IDs to results
    """
    executor = HierarchicalExecutor(num_workers=num_workers, fail_fast=fail_fast)

    # Add all tasks
    for task_def in tasks_definitions:
        executor.add_task(
            task_id=task_def['id'],
            func=task_def['func'],
            args=task_def.get('args', ()),
            kwargs=task_def.get('kwargs', {}),
            dependencies=task_def.get('dependencies', [])
        )

    # Execute and return results
    return executor.execute()
