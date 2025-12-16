#!/bin/bash
#
# run_benchmarks.sh - Script to run CFF Builder benchmarks
#
# Usage:
#   ./run_benchmarks.sh              - Run all predefined benchmarks
#   ./run_benchmarks.sh custom       - Run custom benchmarks (edit CUSTOM_TESTS)
#   ./run_benchmarks.sh single <args> - Run a single benchmark with specific arguments
#
# Examples:
#   ./run_benchmarks.sh single p f 2 1
#   ./run_benchmarks.sh single p g 2 4 1 1
#   ./run_benchmarks.sh single m g 1 2 4 1 1

# Configuration
EXECUTABLE="./generate_cff_benchmark"
OUTPUT_FILE="benchmark_results.txt"

# PREDEFINED TESTS - Edit this section to add/remove tests

# Format: "type arguments"
# Examples:
#   "p f 2 1"           -> ./generate_cff benchmark p f 2 1
#   "p g 2 4 1 1"       -> ./generate_cff benchmark p g 2 4 1 1
#   "m g 1 2 4 1 1"     -> ./generate_cff benchmark m g 1 2 4 1 1

PREDEFINED_TESTS=(
    "p f 2 1"
    "p g 2 4 1 1"
    "p g 4 4 1 2"
    "p g 4 16 2 2"
    "p g 16 16 2 3"
    "p g 16 16 3 4"
)

# Custom tests (used with ./run_benchmarks.sh custom)
CUSTOM_TESTS=(
    "p f 3 1"
    "p g 3 9 1 1"
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

run_single_benchmark() {
    local args="$@"
    echo "Running: $EXECUTABLE benchmark $args"
    $EXECUTABLE benchmark $args
}

run_predefined_tests() {
    local test_count=${#PREDEFINED_TESTS[@]}
    local current=1
    
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

run_custom_tests() {
    local test_count=${#CUSTOM_TESTS[@]}
    local current=1
    
    echo "Total tests: $test_count"
    echo ""
    
    for test in "${CUSTOM_TESTS[@]}"; do
        echo "--------------------------------------------------------------------------------"
        echo "Test $current of $test_count"
        echo "--------------------------------------------------------------------------------"
        
        $EXECUTABLE benchmark $test
        
        ((current++))
    done
}

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
    echo "  $0              - Run all predefined benchmarks"
    echo "  $0 custom       - Run custom benchmarks"
    echo "  $0 single <args> - Run a single benchmark"
    echo "  $0 save         - Run and save results to file"
    echo ""
    echo "Examples:"
    echo "  $0 single p f 2 1"
    echo "  $0 single p g 2 4 1 1"
    echo "  $0 single m g 1 2 4 1 1"
}

# MAIN

print_header
check_executable

case "${1:-}" in
    "")
        echo "Running predefined tests..."
        run_predefined_tests
        ;;
    "custom")
        echo "Running custom tests..."
        run_custom_tests
        ;;
    "single")
        shift
        if [ $# -lt 3 ]; then
            echo "Error: Insufficient arguments for single test."
            print_usage
            exit 1
        fi
        run_single_benchmark "$@"
        ;;
    "save")
        echo "Running tests and saving results..."
        run_predefined_tests
        save_results
        ;;
    "help"|"-h"|"--help")
        print_usage
        ;;
    *)
        echo "Unknown option: $1"
        print_usage
        exit 1
        ;;
esac

echo ""
echo "================================================================================"
echo "                            BENCHMARK COMPLETE                                  "
echo "================================================================================"