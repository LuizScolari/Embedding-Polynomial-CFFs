# Nome do executável
TARGET = update_cff

# Compilador (Apple Clang)
CC = gcc

# --- CAMINHOS DO HOMEBREW (Automático) ---
# Detecta onde o brew está instalado (/opt/homebrew no M1/M2/M3)
BREW_PREFIX := $(shell brew --prefix)

# --- CONFIGURAÇÃO DO OPENMP PARA CLANG ---
# O Clang precisa de -Xpreprocessor e do caminho explícito da libomp
OMP_CFLAGS = -Xpreprocessor -fopenmp -I$(BREW_PREFIX)/opt/libomp/include
OMP_LDFLAGS = -L$(BREW_PREFIX)/opt/libomp/lib -lomp

# --- CONFIGURAÇÃO DA GLIB ---
GLIB_CFLAGS = $(shell pkg-config --cflags glib-2.0)
GLIB_LDFLAGS = $(shell pkg-config --libs glib-2.0)

# --- FLAGS DE COMPILAÇÃO ---
# Inclui headers do sistema + libomp + glib
CFLAGS = -g -Wall -Wextra \
         -I$(BREW_PREFIX)/include \
         $(OMP_CFLAGS) \
         $(GLIB_CFLAGS)

# --- FLAGS DE LINKAGEM ---
# Linka libomp + flint + gmp + math + glib
LDFLAGS = -L$(BREW_PREFIX)/lib \
          -lflint -lgmp -lm \
          $(OMP_LDFLAGS) \
          $(GLIB_LDFLAGS)

# Arquivos fonte
SOURCES = main.c cff_builder.c cff_file_generator.c
OBJECTS = $(SOURCES:.c=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)