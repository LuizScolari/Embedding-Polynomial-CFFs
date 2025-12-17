# CFF Builder Makefile
# (Support for Linux and macOS)

TARGET = generate_cff

# Compiler
CC = gcc

# OPERATING SYSTEM DETECTION
UNAME_S := $(shell uname -s)

# GLIB CONFIGURATION
# (Note: On Linux, ensure libglib2.0-dev is installed)
GLIB_CFLAGS = $(shell pkg-config --cflags glib-2.0)
GLIB_LDFLAGS = $(shell pkg-config --libs glib-2.0)

# PLATFORM SPECIFIC CONFIGURATION
ifeq ($(UNAME_S),Darwin)
    # macOS Configuration (Homebrew + Clang)
    BREW_PREFIX := $(shell brew --prefix)
    
    # OpenMP for Apple Clang
    OMP_CFLAGS = -Xpreprocessor -fopenmp -I$(BREW_PREFIX)/opt/libomp/include
    OMP_LDFLAGS = -L$(BREW_PREFIX)/opt/libomp/lib -lomp
    
    # Homebrew Include/Lib paths
    PLATFORM_INCLUDES = -I$(BREW_PREFIX)/include
    PLATFORM_LIBS = -L$(BREW_PREFIX)/lib
else
    # Linux Configuration (Standard GCC)
    # OpenMP for GCC
    OMP_CFLAGS = -fopenmp
    OMP_LDFLAGS = -fopenmp
    
    # Standard paths (usually empty as libs are in /usr/lib)
    PLATFORM_INCLUDES =
    PLATFORM_LIBS =
endif

# COMPILATION FLAGS
CFLAGS = -g -Wall -Wextra \
         $(PLATFORM_INCLUDES) \
         $(OMP_CFLAGS) \
         $(GLIB_CFLAGS)

# LINKING FLAGS
LDFLAGS = $(PLATFORM_LIBS) \
          -lflint -lgmp -lm \
          $(OMP_LDFLAGS) \
          $(GLIB_LDFLAGS)

# Directories
SRC_DIR = src
BUILD_DIR = build

# Source files
SOURCES = $(SRC_DIR)/main.c $(SRC_DIR)/cff_builder.c $(SRC_DIR)/cff_file_generator.c
OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

.PHONY: all clean dirs

all: dirs $(TARGET)

dirs:
	@mkdir -p $(BUILD_DIR)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)