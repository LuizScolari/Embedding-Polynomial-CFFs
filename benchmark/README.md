# CFF Builder - Automated Benchmark System

## Description

This benchmark system automatically measures the execution time of CFF Builder with two metrics:

- **Time 1**: Total time for inverted index + CFF generation (+ concatenation for embedding)
- **Time 2**: Time for CFF matrix generation + concatenation only

Each test is executed **10 times** and the final result is the **average** of the times.

## Compilation

```bash
make
```

## Usage

### Automatic Benchmark Mode (all predefined tests)

```bash
./generate_cff_benchmark benchmark
```

This automatically executes:
1. `./generate_cff p f 2 1` (10 times, calculates times, saves last CFF)
2. `./generate_cff p g 2 4 1 1` (10 times, calculates times, saves last CFF)

### Individual Benchmark

```bash
# Polynomial from scratch
./generate_cff_benchmark benchmark p f <q> <k>

# Polynomial embedding
./generate_cff_benchmark benchmark p g <q0> <q1> <k0> <k1>

# Monotone embedding
./generate_cff_benchmark benchmark m g <d> <q0> <q1> <k0> <k1>
```

**Examples:**
```bash
./generate_cff_benchmark benchmark p f 2 1
./generate_cff_benchmark benchmark p g 2 4 1 1
./generate_cff_benchmark benchmark p g 4 8 1 1
./generate_cff_benchmark benchmark m g 1 2 4 1 1
```

### Normal Mode (without benchmark)

The program still works normally without benchmark:

```bash
./generate_cff_benchmark p f 2 1
./generate_cff_benchmark p g 2 4 1 1
```

### Benchmark Script (optional)

Use the `run_benchmarks.sh` script for more flexibility:

```bash
# Run predefined tests
./run_benchmarks.sh

# Run custom tests (edit CUSTOM_TESTS array in script)
./run_benchmarks.sh custom

# Run a specific test
./run_benchmarks.sh single p f 2 1

# Save results to file
./run_benchmarks.sh save
```

## Test Configuration

### In main.c

To change predefined tests in `benchmark` mode, edit the `run_all_benchmarks()` function in `main.c`:

```c
void run_all_benchmarks(void) {
    /* Test 1: ./generate_cff p f 2 1 */
    run_benchmark_f('p', 2, 1);
    
    /* Test 2: ./generate_cff p g 2 4 1 1 */
    long Fq_steps_test2[2] = {2, 4};
    long k_steps_test2[2] = {1, 1};
    run_benchmark_g('p', 0, Fq_steps_test2, k_steps_test2);
    
    /* Add more tests here... */
}
```

### In run_benchmarks.sh

Edit the `PREDEFINED_TESTS` or `CUSTOM_TESTS` arrays:

```bash
PREDEFINED_TESTS=(
    "p f 2 1"
    "p g 2 4 1 1"
    "p g 4 8 1 1"
)

CUSTOM_TESTS=(
    "p f 3 1"
    "p g 3 9 1 1"
)
```

### Number of Iterations

To change the number of iterations (default: 10), edit `main.c`:

```c
#define BENCHMARK_ITERATIONS 10
```

## Example Output

```
================================================================================
BENCHMARK: ./generate_cff p f 2 1
Running 10 iterations...
================================================================================
  Iteration 1/10...
  Iteration 2/10...
  ...
  Iteration 10/10...
Generating initial CFF in 'CFFs/1-CFF(4,4).txt'...
Initial CFF matrix of 4x4 generated.

--------------------------------------------------------------------------------
RESULTS for: ./generate_cff p f 2 1
--------------------------------------------------------------------------------
  Iterations: 10
  Time 1 (Inverted Index + CFF Generation):  0.000234 seconds (average)
  Time 2 (Only CFF Matrix Generation):       0.000089 seconds (average)
--------------------------------------------------------------------------------
```

## What Each Time Measures

### In `generate_cff` (action 'f'):
- **Time 1**: Complete `generate_new_cff_blocks()`
- **Time 2**: Only `generate_single_cff()` (internal)

### In `embeed_cff` (action 'g'):
- **Time 1**: `generate_new_cff_blocks()` + final matrix allocation + 4 concatenation loops
- **Time 2**: 3 calls to `generate_single_cff()` + 4 concatenation loops

## Files

| File | Description |
|------|-------------|
| `cff_builder.c` | Main code with time instrumentation |
| `cff_builder.h` | Header with benchmark declarations |
| `cff_file_generator.c` | File read/write functions |
| `cff_file_generator.h` | File generator header |
| `main.c` | Entry point with benchmark system |
| `Makefile` | Build script |
| `run_benchmarks.sh` | Auxiliary script for benchmarks |

## Dependencies

- GCC with OpenMP support
- FLINT (Fast Library for Number Theory)
- GMP (GNU Multiple Precision)
- MPFR
- GLib 2.0

## Clean

```bash
make clean
```