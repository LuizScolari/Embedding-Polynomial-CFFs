# Makefile for CFF Builder with Benchmark Support
# (Compatible with Linux and macOS)
# 
# Usage:
#   make              - Build the benchmark version
#   make clean        - Remove compiled files
#   make run          - Run all benchmarks
#   make run-custom   - Run with custom parameters (edit ARGS variable)

TARGET = generate_cff_benchmark

# Compiler
CC = gcc

# OPERATING SYSTEM DETECTION
UNAME_S := $(shell uname -s)

#  GLIB CONFIGURATION
# Note: On Linux, ensure libglib2.0-dev is installed
GLIB_CFLAGS = $(shell pkg-config --cflags glib-2.0)
GLIB_LDFLAGS = $(shell pkg-config --libs glib-2.0)

# PLATFORM SPECIFIC CONFIGURATION 
ifeq ($(UNAME_S),Darwin)
    # macOS Configuration (Homebrew) 
    BREW_PREFIX := $(shell brew --prefix)
    
    # OpenMP for Clang/Apple
    OMP_CFLAGS = -Xpreprocessor -fopenmp -I$(BREW_PREFIX)/opt/libomp/include
    OMP_LDFLAGS = -L$(BREW_PREFIX)/opt/libomp/lib -lomp
    
    # Include/Lib paths for Homebrew
    PLATFORM_INCLUDES = -I$(BREW_PREFIX)/include
    PLATFORM_LIBS = -L$(BREW_PREFIX)/lib
else
    # Linux Configuration (Standard GCC)
    # OpenMP for GCC
    OMP_CFLAGS = -fopenmp
    OMP_LDFLAGS = -fopenmp
    
    # Standard paths (usually empty if libs are in /usr/lib)
    PLATFORM_INCLUDES =
    PLATFORM_LIBS =
endif

# COMPILATION FLAGS
CFLAGS = -O3 -Wall -Wextra \
         $(PLATFORM_INCLUDES) \
         $(OMP_CFLAGS) \
         $(GLIB_CFLAGS)

# LINKING FLAGS
LDFLAGS = $(PLATFORM_LIBS) \
          -lflint -lgmp -lmpfr -lm \
          $(OMP_LDFLAGS) \
          $(GLIB_LDFLAGS)

# Source files
SRCS = main_benchmark.c cff_builder_benchmark.c cff_file_generator.c
HDRS = cff_builder_benchmark.h cff_file_generator.h

# Phony targets
.PHONY: all clean run run-custom

# Default target
all: $(TARGET)
	@echo ""
	@echo "Build complete: $(TARGET)"
	@echo ""
	@echo "Usage:"
	@echo "  ./$(TARGET) benchmark                           - Run all predefined benchmarks"
	@echo "  ./$(TARGET) benchmark p f <q> <k>               - Single benchmark from scratch"
	@echo "  ./$(TARGET) benchmark p g <q0> <q1> <k0> <k1>   - Single benchmark embedding"
	@echo "  ./$(TARGET) p f <q> <k>                         - Normal mode (no benchmark)"
	@echo "  ./$(TARGET) p g <q0> <q1> <k0> <k1>             - Normal mode embedding"

$(TARGET): $(SRCS) $(HDRS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

# Run all benchmarks
run: $(TARGET)
	./$(TARGET) benchmark

# Run with custom arguments (edit ARGS as needed)
ARGS = benchmark p f 2 1
run-custom: $(TARGET)
	./$(TARGET) $(ARGS)

# Clean build artifacts
clean:
	rm -f $(TARGET)
	rm -rf CFFs/