# CFF Builder - Automated Benchmark System

## Description

This benchmark system automatically measures the execution time of CFF Builder with two metrics:

- **Time 1**: Total time for inverted index + CFF generation (+ concatenation for embedding)
- **Time 2**: Time for CFF matrix generation + concatenation only

Each test is executed **100 times** by default and the final result is the **average** of the times.

> **Important:** To generate a larger CFF, the previous CFF must be saved in a file. For this reason, the tests are executed in construction sequence order. However, file reading and writing operations are **not** included in the measured time — only the CFF generation algorithm is timed.

> **Note:** The benchmark tests are computationally intensive. To reduce execution time, you can change the number of iterations in `main_benchmark.c` by editing the `BENCHMARK_ITERATIONS` constant.

## Compilation

```bash
make
```

## Usage

### Quick Start (Predefined Fast Tests)

```bash
./run_benchmarks.sh
```

This runs a small set of fast tests for quick validation.

### Paper Benchmark Tests

The following test suites were used in the paper and can be run with the benchmark script:

| Command | Description | Complexity |
|---------|-------------|------------|
| `./run_benchmarks.sh f2` | Tests over F₂ field | Heavy |
| `./run_benchmarks.sh f3` | Tests over F₃ field | Heavy |
| `./run_benchmarks.sh f5` | Tests over F₅ field | Moderate |
| `./run_benchmarks.sh monotone` | Monotone CFF tests | Moderate |

**Warning:** These tests are computationally expensive and may take a long time to complete.

### Individual Benchmark

Using the script:
```bash
./run_benchmarks.sh single p f 2 1
./run_benchmarks.sh single p g 2 4 1 1
./run_benchmarks.sh single m g 1 2 4 1 1
```

Or directly with the executable:
```bash
# Initial CFF
./generate_cff_benchmark benchmark p f <q> <k>

# Embedding CFFs
./generate_cff_benchmark benchmark p g <q0> <q1> <k0> <k1>

# Monotone CFFs
./generate_cff_benchmark benchmark m g <d> <q0> <q1> <k0> <k1>
```

### Save Results to File

```bash
./run_benchmarks.sh save
```

Runs the predefined fast tests and saves results to `benchmark_results.txt`.

### Normal Mode (without benchmark)

The program also works without benchmark (single execution):

```bash
./generate_cff_benchmark p f 2 1
./generate_cff_benchmark p g 2 4 1 1
./generate_cff_benchmark m g 1 2 4 1 1
```

## Test Suites Detail

### F₂ Field Tests

```bash
./run_benchmarks.sh f2
```

Tests: 

`p f 2 1` → `p g 2 4 1 1` → `p g 4 4 1 2` → `p g 4 16 2 2` → `p g 16 16 2 3` → `p g 16 16 3 4`

### F₃ Field Tests

```bash
./run_benchmarks.sh f3
```

Tests: 

`p f 3 1` → `p g 3 3 1 2` → `p g 3 9 2 2` → `p g 9 9 2 3` → `p g 9 9 3 4` → `p g 9 9 4 5`

### F₅ Field Tests

```bash
./run_benchmarks.sh f5
```

Tests: 

`p f 5 1` → `p g 5 5 1 2` → `p g 5 25 2 2` → `p g 25 25 2 3`

### Monotone Tests

```bash
./run_benchmarks.sh monotone
```

Tests: 

`p f 2 1` → `m g 1 2 4 1 1` → `m g 1 4 16 1 1` → `m g 1 16 256 1 1`  

`p f 3 1` → `m g 2 3 9 1 1` → `m g 2 9 81 1 1` 

`p f 3 2` → `m g 1 3 27 2 2` 

`p f 5 2` → `m g 2 5 25 2 2`

## Configuration

### Changing Number of Iterations

To change the number of iterations (default: 100), edit `main_benchmark.c`:

```c
#define BENCHMARK_ITERATIONS 100  // Change this value
```

## Example Output

```
================================================================================
  TEST: ./generate_cff p f 2 1
  Running 100 iterations...
================================================================================

  RESULTS:
    Iterations: 100
    Time 1 (Inverted Index + CFF Generation):  0.000147 seconds (average)
    Time 2 (Only CFF Matrix Generation):       0.000043 seconds (average)
================================================================================

================================================================================
  TEST: ./generate_cff p g 2 4 1 1
  Running 100 iterations...
================================================================================

  RESULTS:
    Iterations: 100
    Time 1 (Inverted Index + CFF Generation + Concatenation): 0.000169 seconds (average)
    Time 2 (Only CFF Matrix Generation + Concatenation):      0.000060 seconds (average)
================================================================================
```

## What Each Time Measures

### In `generate_cff` (action 'f'):
- **Time 1**: Complete `generate_new_cff_blocks()`
- **Time 2**: Only `generate_single_cff()` (internal)

### In `embeed_cff` (action 'g'):
- **Time 1**: `generate_new_cff_blocks()` + 4 concatenation loops
- **Time 2**: 3 calls to `generate_single_cff()` + 4 concatenation loops

## Files

| File | Description |
|------|-------------|
| `cff_builder_benchmark.c` | Main code with time instrumentation |
| `cff_builder_benchmark.h` | Header with benchmark declarations |
| `cff_file_generator.c` | File read/write functions |
| `cff_file_generator.h` | File generator header |
| `main_benchmark.c` | Entry point with benchmark system |
| `Makefile` | Build script |
| `run_benchmarks.sh` | Benchmark automation script |

## Dependencies

- GCC with OpenMP support
- FLINT 3.3.1 (Fast Library for Number Theory)
- GMP (GNU Multiple Precision)
- MPFR
- GLib 2.0

## Clean

```bash
make clean
```
