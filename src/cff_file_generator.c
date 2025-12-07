#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <stdint.h>
#include "cff_file_generator.h"
#include "cff_builder.h"

int* parse_int_list(char* str, int* count);

struct cff_parameters* read_parameters(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Arquivo de entrada '%s' não encontrado.\n", filename);
        return NULL;
    }

    struct cff_parameters* params = (struct cff_parameters*) malloc(sizeof(struct cff_parameters));
    if (params == NULL) {
        printf("Erro: Falha ao alocar memória para parâmetros.\n");
        fclose(file);
        return NULL;
    }

    params->Fqs = NULL;
    params->ks = NULL;
    params->fqs_count = 0;
    params->ks_count = 0;

    char* line = NULL;
    size_t len = 0;

    if (getline(&line, &len, file) == -1) {
        printf("Erro: Arquivo está vazio ou não foi possível ler a linha.\n");
        free(params);
        fclose(file);
        if (line) free(line);
        return NULL;
    }

    char fqs_str[256]; 
    char ks_str[256];

    int items_scanned = sscanf(line, "%c [%[^]]] [%[^]]]", 
                               &params->construction, fqs_str, ks_str);

    if (items_scanned != 3) {
        printf("Erro: Formato da primeira linha inválido.\n");
        free(params);
        fclose(file);
        free(line);
        return NULL;
    }

    params->Fqs = parse_int_list(fqs_str, &params->fqs_count);
    params->ks  = parse_int_list(ks_str,  &params->ks_count);

    fclose(file);
    free(line);

    if ((params->fqs_count > 0 && params->Fqs == NULL) || 
        (params->ks_count > 0 && params->ks == NULL)) {
        printf("Erro: Falha ao alocar memória para as listas Fqs/ks.\n");
        free(params->Fqs); 
        free(params->ks);
        free(params);
        return NULL;
    }
    
    return params;
}

uint64_t** read_cff_from_file(const char* filename, long* rows, long* cols) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Arquivo de entrada '%s' não encontrado. Começando do zero.\n", filename);
        *rows = 0; *cols = 0;
        return NULL;
    }

    *rows = 0; *cols = 0;
    char* line = NULL;
    size_t len = 0;

    //Pula a primeira linha
    if (getline(&line, &len, file) == -1) {
        fclose(file);
        if(line) free(line);
        printf("Arquivo '%s' está vazio ou não tem dados após a primeira linha.\n", filename);
        return NULL;
    }
    
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
    getline(&line, &len, file);
    
    long words_per_row = WORDS_FOR_BITS(*cols);
    uint64_t** matrix = (uint64_t**) malloc(*rows * sizeof(uint64_t*));
    
    for (long i = 0; i < *rows; i++) {
        matrix[i] = (uint64_t*) calloc(words_per_row, sizeof(uint64_t));
        for (long j = 0; j < *cols; j++) {
            int bit;
            if (fscanf(file, "%d", &bit) == 1 && bit == 1) {
                SET_BIT(matrix[i], j);
            }
        }
    }

    fclose(file);
    if(line) free(line);
    return matrix;
}

void write_cff_to_file(const char* filename, char construction, long* Fq_steps, int fqs_count, long* K_steps, int ks_count,  uint64_t** matrix, long rows, long cols){
    FILE* file = fopen(filename, "w"); 
    if (file == NULL) {
        printf("Erro ao abrir o arquivo '%s' para escrita.\n", filename);
        return;
    }

    fprintf(file, "%c ", construction);

    fprintf(file, "[");
    for (int i = 0; i < fqs_count; i++) {
        fprintf(file, "%ld", Fq_steps[i]);
        if (i < fqs_count - 1) {
            fprintf(file, ",");
        }
    }
    fprintf(file, "] "); 

    fprintf(file, "[");
    for (int i = 0; i < ks_count; i++) {
        fprintf(file, "%ld", K_steps[i]);
        if (i < ks_count - 1) {
            fprintf(file, ",");
        }
    }
    fprintf(file, "]\n"); 

    for (long i = 0; i < rows; i++) {
        for (long j = 0; j < cols; j++) {
            fprintf(file, "%d", GET_BIT(matrix[i], j) ? 1 : 0);
            if (j < cols - 1) {
                fprintf(file, " "); 
            }
        }
        fprintf(file, "\n"); 
    }

    fclose(file);
}

// Função auxiliar para parse dos parametros
int* parse_int_list(char* str, int* count) {
    *count = 0;
    
    if (str == NULL || strlen(str) == 0) {
        return NULL;
    }

    int commas = 0;
    for (int i = 0; str[i] != '\0'; i++) {
        if (str[i] == ',') {
            commas++;
        }
    }
    *count = commas + 1;

    int* list = (int*) malloc(*count * sizeof(int));
    if (list == NULL) {
        *count = 0;
        return NULL; 
    }

    int i = 0;
    char* token = strtok(str, ",");
    
    while (token != NULL) {
        list[i] = atoi(token);
        i++;
        token = strtok(NULL, ",");
    }

    return list;
}