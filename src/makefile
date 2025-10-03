# Nome do programa executável final
TARGET = update_cff

# Compilador
CC = gcc

# --- Flags da GLib obtidas via pkg-config ---
GLIB_CFLAGS = $(shell pkg-config --cflags glib-2.0)
GLIB_LDFLAGS = $(shell pkg-config --libs glib-2.0)

# Flags de compilação: -g (debug), -Wall (todos os avisos), -Wextra (mais avisos), -I (onde achar os headers)
CFLAGS = -g -Wall -Wextra -I/opt/homebrew/include $(GLIB_CFLAGS)

# Flags de linkagem: -L (onde achar as bibliotecas), -l (quais bibliotecas linkar)
LDFLAGS = -L/opt/homebrew/lib -lflint -lgmp -lm $(GLIB_LDFLAGS)

# Arquivos fonte e objeto
SOURCES = main.c cff_builder.c cff_file_generator.c
OBJECTS = $(SOURCES:.c=.o)

# --- .PHONY para alvos que não são arquivos ---
.PHONY: all clean

# Regra principal: compila tudo
all: $(TARGET)

# Regra de linkagem: como criar o executável a partir dos arquivos objeto
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

# Regra de compilação: como criar um arquivo .o a partir de um .c
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Regra de limpeza: remove os arquivos gerados
clean:
	rm -f $(OBJECTS) $(TARGET)