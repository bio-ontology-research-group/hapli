"""
Hierarchical executor module.

Implements parallel execution that respects hierarchical dependencies
using a directed acyclic graph to manage execution order.
"""

import logging
import time
import networkx as nx
from typing import Dict, List, Set, Callable, Any, Optional, Tuple, Union
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed

from src.parallel.task_manager import create_worker_pool

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
    error: Exception = None
    
    def __post_init__(self):
        """Initialize default values."""
        if self.kwargs is None:
            self.kwargs = {}
        if self.dependencies is None:
            self.dependencies = []
    
    def execute(self) -> Any:
        """
        Execute the task.
        
        Returns:
            Result of the task
            
        Raises:
            Any exception raised during task execution
        """
        self.status = 'running'
        try:
            self.result = self.func(*self.args, **(self.kwargs or {}))
            self.status = 'completed'
            return self.result
        except Exception as e:
            self.status = 'failed'
            self.error = e
            raise e


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
        self.tasks = {}  # Dict[str, Task]
        self.dag = nx.DiGraph()
        self.results = {}  # Dict[str, Any]
        self.errors = {}  # Dict[str, Exception]
        
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
            if dep_id not in self.tasks:
                raise ValueError(f"Unknown dependency: {dep_id} for task {task_id}")
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
        # Check for cycles
        if not nx.is_directed_acyclic_graph(self.dag):
            raise nx.NetworkXUnfeasible("Dependency graph contains cycles")
            
        # Clear previous results
        self.results = {}
        self.errors = {}
        
        # Get topological order (ensures dependencies are executed first)
        execution_order = list(nx.topological_sort(self.dag))
        
        # Set up worker pool
        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            # Keep track of running futures
            futures = {}
            
            # Process tasks in topological order
            while execution_order or futures:
                # Submit tasks whose dependencies are satisfied
                self._submit_ready_tasks(execution_order, futures, executor)
                
                # Process completed futures
                if not self._process_completed_futures(futures):
                    break  # Stop if fail_fast and there was an error
                    
                # Small sleep to prevent CPU hammering in the loop
                if futures:
                    time.sleep(0.01)
                    
        return self.results
    
    def _submit_ready_tasks(self, 
                           remaining_tasks: List[str],
                           futures: Dict[str, Any],
                           executor: ThreadPoolExecutor) -> None:
        """
        Submit tasks whose dependencies are satisfied.
        
        Args:
            remaining_tasks: List of tasks not yet submitted
            futures: Dictionary mapping task IDs to their futures
            executor: ThreadPoolExecutor instance
        """
        i = 0
        while i < len(remaining_tasks):
            task_id = remaining_tasks[i]
            task = self.tasks[task_id]
            
            # Check if all dependencies are completed
            deps_satisfied = all(
                dep_id in self.results and dep_id not in self.errors
                for dep_id in task.dependencies
            )
            
            if deps_satisfied:
                # Remove from remaining tasks
                remaining_tasks.pop(i)
                
                # Instead of trying to inspect functions (which can fail with lambdas),
                # we'll use a simpler approach where we prepare all dependency results
                # but only pass them if they're actually used by the function
                dependency_results = {
                    f"dep_{dep_id}": self.results[dep_id] 
                    for dep_id in task.dependencies
                }
                
                # Create a function that includes dependency results in its closure
                def execute_with_deps():
                    # Copy the task's kwargs to avoid modifying the original
                    task_kwargs = task.kwargs.copy() if task.kwargs else {}
                    
                    # Try to execute with dependency results added to kwargs
                    try:
                        return task.func(*task.args, **{**task_kwargs, **dependency_results})
                    except TypeError:
                        # If that fails (likely due to unexpected kwargs), 
                        # try with just the original kwargs
                        return task.func(*task.args, **task_kwargs)
                
                # Submit the wrapper function that handles dependencies appropriately
                future = executor.submit(execute_with_deps)
                futures[task_id] = future
            else:
                i += 1
    
    def _process_completed_futures(self, futures: Dict[str, Any]) -> bool:
        """
        Process completed futures and update results.
        
        Args:
            futures: Dictionary mapping task IDs to their futures
            
        Returns:
            False if execution should stop (fail_fast and error occurred), True otherwise
        """
        # Check for completed futures
        for task_id, future in list(futures.items()):
            if future.done():
                # Remove from futures dict
                del futures[task_id]
                task = self.tasks[task_id]
                
                try:
                    # Get result and store
                    result = future.result()
                    self.results[task_id] = result
                    logger.debug(f"Task completed: {task_id}")
                except Exception as e:
                    # Record error
                    self.errors[task_id] = e
                    logger.error(f"Task failed: {task_id} - {str(e)}")
                    
                    # Handle fail_fast behavior
                    if self.fail_fast:
                        logger.error(f"Stopping execution due to task failure: {task_id}")
                        for f in futures.values():
                            f.cancel()
                        return False
                        
        return True
    
    def get_execution_plan(self) -> List[List[str]]:
        """
        Get the execution plan showing which tasks can execute in parallel.
        
        Returns:
            List of lists, where each inner list contains tasks that can
            execute in parallel at the same level.
        """
        # Check for cycles
        if not nx.is_directed_acyclic_graph(self.dag):
            raise nx.NetworkXUnfeasible("Dependency graph contains cycles")
            
        # Find longest path from any root to each node
        # This determines the "level" of each task
        levels = {}
        for node in self.dag.nodes():
            # If node has no predecessors, it's a root at level 0
            if not list(self.dag.predecessors(node)):
                levels[node] = 0
                
        # Do a topological traversal to find levels
        for node in nx.topological_sort(self.dag):
            node_level = levels.get(node, 0)
            for succ in self.dag.successors(node):
                levels[succ] = max(levels.get(succ, 0), node_level + 1)
                
        # Group tasks by level
        max_level = max(levels.values()) if levels else 0
        execution_plan = [[] for _ in range(max_level + 1)]
        
        for node, level in levels.items():
            execution_plan[level].append(node)
            
        return execution_plan


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
