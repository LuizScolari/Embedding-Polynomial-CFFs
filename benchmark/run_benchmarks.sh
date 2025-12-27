#!/bin/bash
#
# run_benchmarks.sh - CFF Builder Benchmark Automation Script
#
# This script automates the execution of CFF Builder benchmarks.
# Each test runs 100 iterations by default (configurable in main_benchmark.c).
#
# USAGE
#
#   ./run_benchmarks.sh              - Run predefined fast tests (quick validation)
#   ./run_benchmarks.sh f2           - Run F2 field tests (paper benchmark)
#   ./run_benchmarks.sh f3           - Run F3 field tests (paper benchmark)
#   ./run_benchmarks.sh f5           - Run F5 field tests (paper benchmark)
#   ./run_benchmarks.sh monotone     - Run monotone tests (paper benchmark)
#   ./run_benchmarks.sh single <args> - Run a single benchmark with specific args
#   ./run_benchmarks.sh save         - Run predefined tests and save to file
#
# EXAMPLES
#
#   ./run_benchmarks.sh single p f 2 1
#   ./run_benchmarks.sh single p g 2 4 1 1
#   ./run_benchmarks.sh single m g 1 2 4 1 1
#
# WARNING
#
# The F2 and F3 test suites are computationally intensive! F5 and MONOTONE are moderate.
# They were used for the paper benchmarks and may take a long time to complete.
# Use the predefined fast tests (no arguments) for quick validation.
#

# CONFIGURATION

EXECUTABLE="./generate_cff_benchmark"
OUTPUT_FILE="benchmark_results.txt"

# TEST SUITES

# PREDEFINED_TESTS: Fast tests for quick validation
# These run quickly and are useful for checking if everything works
PREDEFINED_TESTS=(
    "p f 2 1"
    "p g 2 4 1 1"
    "p g 4 4 1 2"
)

# PAPER BENCHMARK TESTS

# F2_TESTS: Benchmarks over the F2 (binary) field
F2_TESTS=(
    "p f 2 1"
    "p g 2 4 1 1"
    "p g 4 4 1 2"
    "p g 4 16 2 2"
    "p g 16 16 2 3"
    "p g 16 16 3 4"
)

# F3_TESTS: Benchmarks over the F3 field
F3_TESTS=(
    "p f 3 1"
    "p g 3 3 1 2"
    "p g 3 9 2 2"
    "p g 9 9 2 3"
    "p g 9 9 3 4"
    "p g 9 9 4 5"
)

# F5_TESTS: Benchmarks over the F5 field
F5_TESTS=(
    "p f 5 1"
    "p g 5 5 1 2"
    "p g 5 25 2 2"
    "p g 25 25 2 3"
)

# MONOTONE: Monotone CFF benchmarks
MONOTONE=(
    "p f 2 1"
    "m g 1 2 4 1 1"
    "m g 1 4 16 1 1"
    "m g 1 16 256 1 1"
    "p f 3 1"
    "m g 2 3 9 1 1"
    "m g 2 9 81 1 1"
    "p f 3 2"
    "m g 1 3 27 2 2"
    "p f 5 2"
    "m g 2 5 25 2 2"
)

# FUNCTIONS

print_header() {
    echo ""
    echo "================================================================================"
    echo "                     CFF BUILDER - BENCHMARK AUTOMATION                         "
    echo "================================================================================"
    echo ""
}

check_executable() {
    if [ ! -f "$EXECUTABLE" ]; then
        echo "Error: Executable '$EXECUTABLE' not found."
        echo "Run 'make' first to compile the project."
        exit 1
    fi
}

# Run a single benchmark with given arguments
run_single_benchmark() {
    local args="$@"
    echo "Running: $EXECUTABLE benchmark $args"
    $EXECUTABLE benchmark $args
}

# Run the predefined fast tests (quick validation)
run_predefined_tests() {
    local test_count=${#PREDEFINED_TESTS[@]}
    local current=1
    
    echo "Running PREDEFINED tests (fast validation)..."
    echo "Total tests: $test_count"
    echo ""
    
    for test in "${PREDEFINED_TESTS[@]}"; do
        echo "--------------------------------------------------------------------------------"
        echo "Test $current of $test_count"
        echo "--------------------------------------------------------------------------------"
        
        $EXECUTABLE benchmark $test
        
        ((current++))
    done
}

# Run F2 field tests (paper benchmark - heavy)
run_F2_tests() {
    local test_count=${#F2_TESTS[@]}
    local current=1
    
    echo "Running F2 FIELD tests (paper benchmark)..."
    echo "WARNING: These tests are computationally heavy!"
    echo "Total tests: $test_count"
    echo ""
    
    for test in "${F2_TESTS[@]}"; do
        echo "--------------------------------------------------------------------------------"
        echo "Test $current of $test_count"
        echo "--------------------------------------------------------------------------------"
        
        $EXECUTABLE benchmark $test
        
        ((current++))
    done
}

# Run F3 field tests (paper benchmark - heavy)
run_F3_tests() {
    local test_count=${#F3_TESTS[@]}
    local current=1
    
    echo "Running F3 FIELD tests (paper benchmark)..."
    echo "WARNING: These tests are computationally heavy!"
    echo "Total tests: $test_count"
    echo ""
    
    for test in "${F3_TESTS[@]}"; do
        echo "--------------------------------------------------------------------------------"
        echo "Test $current of $test_count"
        echo "--------------------------------------------------------------------------------"
        
        $EXECUTABLE benchmark $test
        
        ((current++))
    done
}

# Run F5 field tests (paper benchmark - moderate)
run_F5_tests() {
    local test_count=${#F5_TESTS[@]}
    local current=1
    
    echo "Running F5 FIELD tests (paper benchmark)..."
    echo "Total tests: $test_count"
    echo ""
    
    for test in "${F5_TESTS[@]}"; do
        echo "--------------------------------------------------------------------------------"
        echo "Test $current of $test_count"
        echo "--------------------------------------------------------------------------------"
        
        $EXECUTABLE benchmark $test
        
        ((current++))
    done
}

# Run monotone tests (paper benchmark - moderate)
run_monotone_tests() {
    local test_count=${#MONOTONE[@]}
    local current=1
    
    echo "Running MONOTONE tests (paper benchmark)..."
    echo "Total tests: $test_count"
    echo ""
    
    for test in "${MONOTONE[@]}"; do
        echo "--------------------------------------------------------------------------------"
        echo "Test $current of $test_count"
        echo "--------------------------------------------------------------------------------"
        
        $EXECUTABLE benchmark $test
        
        ((current++))
    done
}

# Run predefined tests and save results to file
save_results() {
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    
    echo ""
    echo "Saving results to '$OUTPUT_FILE'..."
    
    {
        echo "================================================================================"
        echo "CFF BUILDER BENCHMARK RESULTS"
        echo "Date: $timestamp"
        echo "================================================================================"
        echo ""
        
        for test in "${PREDEFINED_TESTS[@]}"; do
            $EXECUTABLE benchmark $test 2>&1
            echo ""
        done
        
    } >> "$OUTPUT_FILE"
    
    echo "Results saved successfully!"
}

print_usage() {
    echo "Usage:"
    echo "  $0                - Run predefined fast tests (quick validation)"
    echo "  $0 f2             - Run F2 field tests (paper benchmark - heavy)"
    echo "  $0 f3             - Run F3 field tests (paper benchmark - heavy)"
    echo "  $0 f5             - Run F5 field tests (paper benchmark - moderate)"
    echo "  $0 monotone       - Run monotone tests (paper benchmark - moderate)"
    echo "  $0 single <args>  - Run a single benchmark with specific arguments"
    echo "  $0 save           - Run predefined tests and save results to file"
    echo ""
    echo "Examples:"
    echo "  $0 single p f 2 1"
    echo "  $0 single p g 2 4 1 1"
    echo "  $0 single m g 1 2 4 1 1"
    echo ""
    echo "Note: F2, F3, F5 and MONOTONE tests are computationally intensive."
    echo "      Each test runs 100 iterations (configurable in main_benchmark.c)."
}

# =============================================================================
# MAIN
# =============================================================================

print_header
check_executable

case "${1:-}" in
    "")
        run_predefined_tests
        ;;
    "f2")
        run_F2_tests
        ;;
    "f3")
        run_F3_tests
        ;;
    "f5")
        run_F5_tests
        ;;
    "monotone")
        run_monotone_tests
        ;;
    "single")
        shift
        if [ $# -lt 3 ]; then
            echo "Error: Insufficient arguments for single test."
            echo ""
            print_usage
            exit 1
        fi
        run_single_benchmark "$@"
        ;;
    "save")
        echo "Running predefined tests and saving results..."
        run_predefined_tests
        save_results
        ;;
    "help"|"-h"|"--help")
        print_usage
        ;;
    *)
        echo "Unknown option: $1"
        echo ""
        print_usage
        exit 1
        ;;
esac

echo ""
echo "================================================================================"
echo "                            BENCHMARK COMPLETE                                  "
echo "================================================================================"