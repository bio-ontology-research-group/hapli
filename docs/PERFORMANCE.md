# Performance Tuning Guidelines

This document provides guidance on tuning the parallelization parameters for optimal performance.

## Choosing Between Processes and Threads

The parallelization module supports both process-based and thread-based parallelization:

- **Process-based parallelization** (`pool_type='process'`):
  - Best for CPU-bound tasks
  - Bypasses Python's Global Interpreter Lock (GIL)
  - Higher memory overhead due to process isolation
  - Better for tasks that don't share a lot of data

- **Thread-based parallelization** (`pool_type='thread'`):
  - Best for I/O-bound tasks
  - Lower memory overhead as memory is shared
  - Limited by Python's GIL for CPU-bound tasks
  - Better for tasks that share large data structures

General guideline: Use processes for computation-heavy tasks and threads for I/O-bound tasks.

## Number of Workers

The optimal number of workers depends on your hardware and the nature of your tasks:

- For **CPU-bound tasks**, a good starting point is the number of CPU cores:
  ```python
  num_workers = os.cpu_count()
  ```

- For **I/O-bound tasks**, you can often use more workers than you have CPU cores:
  ```python
  num_workers = os.cpu_count() * 2  # or higher
  ```

- If your tasks use significant memory, you may need to reduce the number of workers to avoid memory pressure:
  ```python
  num_workers = max(1, os.cpu_count() // 2)
  ```

## Task Chunking

Proper task chunking is crucial for efficient parallelization:

- **Too small chunks**: High overhead from task distribution and result collection
- **Too large chunks**: Poor load balancing and reduced parallelism

Guidelines for chunk size:

1. For many small tasks, use larger chunks:
   ```python
   chunk_size = max(1, len(tasks) // (num_workers * 10))
   ```

2. For fewer large tasks, use smaller chunks:
   ```python
   chunk_size = max(1, len(tasks) // (num_workers * 2))
   ```

3. A general heuristic:
   ```python
   chunk_size = max(1, min(100, len(tasks) // num_workers))
   ```

## Batch Size for Alignment Tasks

When using `ParallelAligner`, the `batch_size` parameter controls how many features are processed in each batch:

- Smaller batch sizes improve load balancing but increase overhead
- Larger batch sizes reduce overhead but may lead to uneven work distribution

Good starting points:
- For small features (<1kb): 50-100 features per batch
- For medium features (1-10kb): 20-50 features per batch
- For large features (>10kb): 5-20 features per batch

## Memory Considerations

Parallel processing can significantly increase memory usage:

1. **Process-based parallelism** duplicates data across processes
2. **Large reference sequences** can multiply memory usage
3. **Result aggregation** can temporarily double memory needs

Strategies to manage memory:
- Process data in batches
- Use memory profiling to identify bottlenecks
- Consider reducing `num_workers` for memory-intensive operations

## Hierarchical Task Execution

For complex task graphs with dependencies:

1. **Task granularity**: Balance between too fine-grained (high overhead) and too coarse-grained (limited parallelism)
2. **Critical path optimization**: Identify and optimize the longest path in your task graph
3. **Resource allocation**: Allocate more resources to tasks on the critical path

## Profiling and Benchmarking

Always measure performance to guide optimization:

```python
import time

# Serial execution
start_time = time.time()
serial_results = [process_item(item) for item in items]
serial_time = time.time() - start_time

# Parallel execution
start_time = time.time()
parallel_results = execute_parallel(process_item, items, num_workers=4)
parallel_time = time.time() - start_time

speedup = serial_time / parallel_time
print(f"Speedup: {speedup:.2f}x")
```

## Common Performance Pitfalls

1. **Small workloads**: Parallelization overhead can exceed benefits for small tasks
2. **Shared resource bottlenecks**: Database connections, file I/O, network
3. **GIL contention**: Thread-based parallelism may not help CPU-bound Python code
4. **Load imbalance**: Tasks with highly variable execution times
5. **Data serialization**: Large objects being passed between processes

## Practical Examples

### Example 1: CPU-bound alignment tasks

```python
from src.parallel.parallel_alignment import ParallelAligner

aligner = ParallelAligner(
    num_workers=os.cpu_count(),  # Use all CPU cores
    batch_size=50,               # Moderate batch size
    pool_type='process'          # Process-based for CPU work
)

results = aligner.align_features(features, reference_seq)
```

### Example 2: I/O-bound data loading

```python
from src.parallel.task_manager import execute_parallel

def load_file(filename):
    with open(filename, 'r') as f:
        return f.read()

# Use threads for I/O-bound work
files_content = execute_parallel(
    load_file, 
    file_list,
    num_workers=os.cpu_count() * 2,  # More workers for I/O
    pool_type='thread'              # Thread-based for I/O
)
```

### Example 3: Hierarchical analysis workflow

```python
from src.parallel.hierarchical_executor import HierarchicalExecutor

executor = HierarchicalExecutor(num_workers=os.cpu_count())

# Add tasks with dependencies
executor.add_task('load_data', load_function)
executor.add_task('parse', parse_function, dependencies=['load_data'])
executor.add_task('filter', filter_function, dependencies=['parse'])
executor.add_task('analyze', analyze_function, dependencies=['filter'])
executor.add_task('visualize', visualize_function, dependencies=['analyze'])

# Execute respecting dependencies
results = executor.execute()
```

## Further Reading

- [Python Multiprocessing Documentation](https://docs.python.org/3/library/multiprocessing.html)
- [Python Threading Documentation](https://docs.python.org/3/library/threading.html)
- [Python Concurrent.futures Documentation](https://docs.python.org/3/library/concurrent.futures.html)
