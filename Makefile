# CFF Builder Makefile
TARGET = generate_cff

# Compiler (Apple Clang)
CC = gcc

# --- HOMEBREW PATHS (Automatic) ---
BREW_PREFIX := $(shell brew --prefix)

# --- OPENMP CONFIGURATION FOR CLANG ---
OMP_CFLAGS = -Xpreprocessor -fopenmp -I$(BREW_PREFIX)/opt/libomp/include
OMP_LDFLAGS = -L$(BREW_PREFIX)/opt/libomp/lib -lomp

# --- GLIB CONFIGURATION ---
GLIB_CFLAGS = $(shell pkg-config --cflags glib-2.0)
GLIB_LDFLAGS = $(shell pkg-config --libs glib-2.0)

# --- COMPILATION FLAGS ---
CFLAGS = -g -Wall -Wextra \
         -I$(BREW_PREFIX)/include \
         $(OMP_CFLAGS) \
         $(GLIB_CFLAGS)

# --- LINKING FLAGS ---
LDFLAGS = -L$(BREW_PREFIX)/lib \
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





