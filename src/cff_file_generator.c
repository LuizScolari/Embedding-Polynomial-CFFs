#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int** read_cff_from_file(const char* filename, long* rows, long* cols) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Arquivo de entrada '%s' não encontrado. Começando do zero.\n", filename);
        *rows = 0; *cols = 0;
        return NULL;
    }

    *rows = 0; *cols = 0;
    char* line = NULL;
    size_t len = 0;
    
    // Primeiro passo: descobre as dimensões
    if (getline(&line, &len, file) != -1) {
        (*rows)++;
        char* token = strtok(line, " \t\n");
        while (token != NULL) {
            if (strlen(token) > 0) (*cols)++;
            token = strtok(NULL, " \t\n");
        }
    }
    while (getline(&line, &len, file) != -1) {
        (*rows)++;
    }
    
    // Segundo passo: lê os dados
    rewind(file);
    int** matrix = (int**) malloc(*rows * sizeof(int*));
    for (long i = 0; i < *rows; i++) {
        matrix[i] = (int*) malloc(*cols * sizeof(int));
        for (long j = 0; j < *cols; j++) {
            fscanf(file, "%d", &matrix[i][j]);
        }
    }

    fclose(file);
    if(line) free(line);
    printf("Matriz de %ldx%ld lida do arquivo '%s'.\n", *rows, *cols, filename);
    return matrix;
}

void write_cff_to_file(const char* filename, int** matrix, long rows, long cols) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("Erro ao abrir arquivo para escrita");
        return;
    }
    for (long i = 0; i < rows; i++) {
        for (long j = 0; j < cols; j++) {
            // Imprime o número
            fprintf(file, "%d", matrix[i][j]);
            // Imprime um espaço apenas se NÃO for o último número da linha
            if (j < cols - 1) {
                fprintf(file, " ");
            }
        }
        fprintf(file, "\n"); // Quebra de linha ao final de cada linha da matriz
    }
    fclose(file);
    printf("Nova matriz de %ldx%ld salva no arquivo '%s'.\n", rows, cols, filename);
}