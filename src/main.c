#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cff_builder.h" // Inclui nossa "biblioteca"

// --- Funções Auxiliares para Arquivos e Matrizes ---

// Lê uma matriz de um arquivo de texto
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

// Escreve uma matriz para um arquivo de texto
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

// Função para liberar uma matriz int** simples
void free_matrix(int** matrix, long rows) {
    if (!matrix) return;
    for (long i = 0; i < rows; i++) free(matrix[i]);
    free(matrix);
}

int main() {
    const char* filename = "cff_matrix.txt";
    
    // 1. LÊ A MATRIZ ANTIGA (CFF_old_old) DO ARQUIVO
    long old_rows = 0, old_cols = 0;
    int** cff_old_old = read_cff_from_file(filename, &old_rows, &old_cols);

    // 2. CHAMA O CÓDIGO PARA GERAR OS NOVOS BLOCOS
    // Parâmetros para a nova extensão (ex: ir de F_4 para F_16)
    long Fq_steps[] = {2,4,16};
    int num_steps = 3;
    long max_poly_degree = 1;

    printf("\nGerando 3 novos blocos da CFF...\n");
    generated_cffs new_blocks = generate_new_cff_blocks(Fq_steps, num_steps, max_poly_degree);

    // 3. CONCATENA AS 4 PARTES
    long new_total_rows = old_rows + new_blocks.rows_new_old;
    long new_total_cols = old_cols + new_blocks.cols_old_new;
    
    int** final_cff = (int**) malloc(new_total_rows * sizeof(int*));
    for (long i = 0; i < new_total_rows; i++) {
        final_cff[i] = (int*) calloc(new_total_cols, sizeof(int));
    }

    printf("Concatenando matrizes para formar uma matriz final de %ldx%ld.\n", new_total_rows, new_total_cols);

    // Copia bloco 1: CFF_old_old
    for(long i=0; i < old_rows; i++) memcpy(final_cff[i], cff_old_old[i], old_cols * sizeof(int));

    // Copia bloco 2: CFF_old_new
    for(long i=0; i < new_blocks.rows_old_new; i++) memcpy(&final_cff[i][old_cols], new_blocks.cff_old_new[i], new_blocks.cols_old_new * sizeof(int));
    
    // Copia bloco 3: CFF_new_old
    for(long i=0; i < new_blocks.rows_new_old; i++) memcpy(final_cff[old_rows + i], new_blocks.cff_new_old[i], new_blocks.cols_new_old * sizeof(int));
    
    // Copia bloco 4: CFF_new
    for(long i=0; i < new_blocks.rows_new; i++) memcpy(&final_cff[old_rows + i][old_cols], new_blocks.cff_new[i], new_blocks.cols_new * sizeof(int));
    

    // 4. SALVA A NOVA MATRIZ NO ARQUIVO
    write_cff_to_file(filename, final_cff, new_total_rows, new_total_cols);

    // 5. LIMPEZA TOTAL DA MEMÓRIA
    printf("\nLimpando toda a memória.\n");
    free_matrix(cff_old_old, old_rows);
    free_generated_cffs(&new_blocks);
    free_matrix(final_cff, new_total_rows);

    return 0;
}