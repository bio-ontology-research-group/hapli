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
        task = Task(id=task_id, func=func, args=args, kwargs=kwargs or {}, dependencies=dependencies or [])
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
            logger.warning("Could not check for cycles in the dependency graph.")

        # Clear previous results
        self.results = {}
        self.errors = {}
        
        logger.info(f"Starting execution of {len(self.tasks)} tasks with fail_fast={self.fail_fast}")

        # Get execution plan - groups of tasks that can be executed in parallel
        execution_plan = self.get_execution_plan()
        
        # Set up worker pool
        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            # Execute tasks level by level (ensures dependencies are executed first)
            for level, task_ids in enumerate(execution_plan):
                logger.debug(f"Executing level {level} with {len(task_ids)} tasks")
                
                # Dictionary to store futures for this level
                futures = {}
                
                # Submit all tasks at this level
                for task_id in task_ids:
                    # Skip if any dependency had errors
                    task = self.tasks[task_id]
                    if any(dep in self.errors for dep in task.dependencies):
                        msg = f"Skipping task {task_id} due to failed dependencies"
                        logger.warning(msg)
                        self.errors[task_id] = RuntimeError(msg)
                        continue
                    
                    # Get dependency results in correct order
                    dependency_results = []
                    try:
                        for dep_id in task.dependencies:
                            if dep_id not in self.results:
                                raise ValueError(f"Dependency {dep_id} has no result. This is a bug.")
                            dependency_results.append(self.results[dep_id])
                    except Exception as e:
                        logger.error(f"Error collecting dependencies for {task_id}: {e}")
                        self.errors[task_id] = RuntimeError(f"Failed to collect dependencies: {e}")
                        continue
                    
                    # Submit task with its dependency results
                    logger.debug(f"Submitting task {task_id} for execution")
                    future = executor.submit(task.execute, tuple(dependency_results))
                    futures[task_id] = future
                
                # Wait for all futures to complete
                for task_id, future in futures.items():
                    try:
                        result = future.result()  # Blocks until task is done
                        self.results[task_id] = result
                        logger.debug(f"Task completed: {task_id}")
                    except Exception as e:
                        logger.error(f"Task failed: {task_id} - {type(e).__name__}: {e}")
                        self.errors[task_id] = e
                        
                        if self.fail_fast:
                            logger.error(f"Fail fast enabled. Stopping execution due to task failure: {task_id}")
                            
                            # Mark remaining tasks in this level as failed
                            for remaining_id, remaining_future in futures.items():
                                if remaining_id != task_id and not remaining_future.done():
                                    remaining_future.cancel()
                                    self.errors[remaining_id] = RuntimeError(f"Cancelled due to fail_fast from {task_id}")
                            
                            # Mark all tasks in future levels as failed
                            for future_level in range(level + 1, len(execution_plan)):
                                for future_task_id in execution_plan[future_level]:
                                    self.errors[future_task_id] = RuntimeError(f"Execution stopped due to fail_fast on task {task_id}")
                            
                            return self.results  # Early return
                
                # If fail_fast is False, continue with the next level
        
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
