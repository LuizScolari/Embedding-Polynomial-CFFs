// Run command -> gcc -I/opt/homebrew/include -L/opt/homebrew/lib test_flint.c -o test_flint -lflint -lgmp
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // Para calloc
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_poly.h"
#include "flint/nmod_poly.h"

// --- Estruturas de Dados ---

// Struct para organizar os resultados da partição para um subcorpo
typedef struct {
    long q; // Tamanho do subcorpo (ex: 2, para F_2)

    fq_nmod_t* all_elements; // Corresponde a F_sets[q]
    long count_all;          
    long capacity_all;       

    fq_nmod_t* only_elements; // Corresponde a only[q]
    long count_only;         
    long capacity_only;      
} subfield_partition;


// --- Protótipos de Funções ---
void get_element_by_arithmetic(fq_nmod_t result, ulong i, fq_nmod_ctx_t ctx);
void add_element_to_list(fq_nmod_t** list, long* count, long* capacity, const fq_nmod_t element, const fq_nmod_ctx_t ctx);
subfield_partition* partition_by_subfields(const long* Fq_steps, int num_steps, const fq_nmod_ctx_t ctx);
fq_nmod_poly_t* generate_all_polynomials(long* poly_count, long max_degree, const fq_nmod_ctx_t ctx);
void generate_recursive(fq_nmod_poly_t* poly_list, long* current_index, fq_nmod_poly_t current_poly, long degree, long max_degree, const fq_nmod_t* elements, long num_elements, const fq_nmod_ctx_t ctx);
fq_nmod_t** create_evaluation_matrix(long num_points, long num_polys, const fq_nmod_t* points, const fq_nmod_poly_t* polys, const fq_nmod_ctx_t ctx);
void free_evaluation_matrix(fq_nmod_t** matrix, long num_points, long num_polys, const fq_nmod_ctx_t ctx);


int main() {
    printf("FLINT version: %s\n\n", FLINT_VERSION);

    // --- Parâmetros ---
    fmpz_t pz;
    fmpz_init(pz);
    fmpz_set_ui(pz, 2);
    long D = 2; // Campo principal será F_q^d
    
    // Subcorpos que queremos encontrar dentro de F_q
    long Fq_steps[] = {4}; 
    int num_steps = 1;

    // --- INICIALIZAÇÃO DO CORPO ---

    // Criar o polinômio irredutível
    nmod_poly_t mod_poly;
    nmod_poly_init(mod_poly, 2);
    // Definir o polinômio, um irredutível para F_q
    nmod_poly_set_coeff_ui(mod_poly, 2, 1);
    nmod_poly_set_coeff_ui(mod_poly, 1, 1);
    nmod_poly_set_coeff_ui(mod_poly, 0, 1);

    // Inicializar o contexto usando nosso polinômio
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_modulus(ctx, mod_poly, "a");
    
    printf("Contexto inicializado para F_%ld^%ld com módulo: ", fmpz_get_ui(pz), D);
    nmod_poly_print(mod_poly);
    printf("\n\n");
    nmod_poly_clear(mod_poly); // Já podemos limpar o polinômio
    
    // --- Execução da Lógica de Particionamento ---
    subfield_partition* partitions = partition_by_subfields(Fq_steps, num_steps, ctx);

    // --- Impressão dos Resultados ---
    for (int i = 0; i < num_steps; i++) {
        printf("--- Subcorpo F_%ld ---\n", partitions[i].q);
        printf("  Total de elementos (F_sets): %ld\n", partitions[i].count_all);
        printf("  Elementos 'Only': %ld\n", partitions[i].count_only);
        
        printf("  Conteúdo de 'All': { ");
        for(long j = 0; j < partitions[i].count_all; j++) {
            fq_nmod_print_pretty(partitions[i].all_elements[j], ctx);
            if (j < partitions[i].count_all - 1) printf(", ");
        }
        printf(" }\n");

        printf("  Conteúdo de 'Only': { ");
        for(long j = 0; j < partitions[i].count_only; j++) {
            fq_nmod_print_pretty(partitions[i].only_elements[j], ctx);
            if (j < partitions[i].count_only - 1) printf(", ");
        }
        printf(" }\n\n");
    }

    // --- Geração dos Polinômios ---
    long max_poly_degree = 1;
    long num_polys = 0;
    fq_nmod_poly_t* all_polys = generate_all_polynomials(&num_polys, max_poly_degree, ctx);

    // --- AVALIAÇÃO DE POLINÔMIOS (usando a nova função) ---
    if (all_polys != NULL && partitions != NULL) {
        printf("\n--- Avaliando Polinômios ---\n");
        
        // Seleciona os pontos e polinômios para avaliar
        subfield_partition Fq_partition = partitions[0]; 
        long num_points_to_eval = Fq_partition.count_all;

        // Chama a função para criar a matriz de avaliação
        fq_nmod_t** evals = create_evaluation_matrix(
            num_points_to_eval, 
            num_polys, 
            Fq_partition.all_elements, 
            all_polys, 
            ctx
        );

        // Imprime a matriz de resultados para verificação
        if (evals != NULL) {
            printf("Matriz de Avaliação [ponto][polinômio]:\n");
            for (long i = 0; i < num_polys; i++) {
                printf("  Ponto "); fq_nmod_print_pretty(Fq_partition.all_elements[i], ctx); printf(": [ ");
                for (long j = 0; j < num_polys; j++) {
                    fq_nmod_print_pretty(evals[i][j], ctx);
                    if (j < num_polys - 1) printf(", ");
                }
                printf(" ]\n");
            }

            // Chama a função para liberar a memória da matriz
            free_evaluation_matrix(evals, num_points_to_eval, num_polys, ctx);
        }
    }
    
    
    // --- Impressão dos Resultados ---
    if (all_polys != NULL) {
        printf("--- %ld Polinômios Gerados (Grau <= %ld) ---\n", num_polys, max_poly_degree);
        for (long i = 0; i < num_polys; i++) {
            printf("P_%ld(a) = ", i);
            // fq_nmod_poly_print_pretty usa 'x' como variável por padrão
            fq_nmod_poly_print_pretty(all_polys[i], "x", ctx);
            printf("\n");
        }

        // --- Liberação de Memória (MUITO IMPORTANTE) ---
        printf("\nLimpando a memória dos polinômios...\n");
        for (long i = 0; i < num_polys; i++) {
            fq_nmod_poly_clear(all_polys[i], ctx);
        }
        free(all_polys);
        printf("Memória liberada.\n");
    }
    

    // --- Liberação de Memória (MUITO IMPORTANTE) ---
    printf("Limpando a memória...\n");
    for (int i = 0; i < num_steps; i++) {
        for (long j = 0; j < partitions[i].capacity_all; j++) fq_nmod_clear(partitions[i].all_elements[j], ctx);
        for (long j = 0; j < partitions[i].capacity_only; j++) fq_nmod_clear(partitions[i].only_elements[j], ctx);
        free(partitions[i].all_elements);
        free(partitions[i].only_elements);
    }
    free(partitions);
    fmpz_clear(pz);
    fq_nmod_ctx_clear(ctx);
    printf("Memória liberada.\n");


    
    return 0;
}


/**
 * @brief Itera por todos os elementos do corpo principal e os particiona
 * de acordo com sua pertinência aos subcorpos especificados.
 */
subfield_partition* partition_by_subfields(const long* Fq_steps, int num_steps, const fq_nmod_ctx_t ctx) {
    fmpz_t order_z;
    fmpz_init(order_z);
    fq_nmod_ctx_order(order_z, ctx);
    long field_size = fmpz_get_si(order_z);
    fmpz_clear(order_z);

    subfield_partition* partitions = malloc(num_steps * sizeof(subfield_partition));
    char* is_seen = calloc(field_size, sizeof(char));

    fq_nmod_t elem, elem_pow_q;
    fmpz_t q_fmpz;
    fq_nmod_init(elem, ctx);
    fq_nmod_init(elem_pow_q, ctx);
    fmpz_init(q_fmpz);

    for (int i = 0; i < num_steps; i++) {
        long q = Fq_steps[i];
        
        partitions[i] = (subfield_partition){.q = q, .count_all = 0, .capacity_all = 0, .count_only = 0, .capacity_only = 0};
        fmpz_set_ui(q_fmpz, q);

        for (long j = 0; j < field_size; j++) {
            get_element_by_arithmetic(elem, j, ctx);

            fq_nmod_pow(elem_pow_q, elem, q_fmpz, ctx);

            if (fq_nmod_equal(elem_pow_q, elem, ctx)) {
                add_element_to_list(&partitions[i].all_elements, &partitions[i].count_all, &partitions[i].capacity_all, elem, ctx);
                
                if (is_seen[j] == 0) {
                    add_element_to_list(&partitions[i].only_elements, &partitions[i].count_only, &partitions[i].capacity_only, elem, ctx);
                    is_seen[j] = 1;
                }
            }
        }
    }

    free(is_seen);
    fq_nmod_clear(elem, ctx);
    fq_nmod_clear(elem_pow_q, ctx);
    fmpz_clear(q_fmpz);

    return partitions;
}

/**
 * @brief Função auxiliar para gerenciar arrays dinâmicos de elementos fq_nmod_t.
 */
void add_element_to_list(fq_nmod_t** list, long* count, long* capacity, const fq_nmod_t element, const fq_nmod_ctx_t ctx) {
    if (*count >= *capacity) {
        *capacity = (*capacity == 0) ? 8 : (*capacity) * 2;
        *list = realloc(*list, (*capacity) * sizeof(fq_nmod_t));
        // IMPORTANTE: Inicializa os novos espaços alocados no array
        for (long i = *count; i < *capacity; i++) {
            fq_nmod_init((*list)[i], ctx);
        }
    }
    fq_nmod_set((*list)[*count], element, ctx);
    (*count)++;
}

/**
 * @brief Constrói o elemento via aritmética.
 */
void get_element_by_arithmetic(fq_nmod_t result, ulong i, fq_nmod_ctx_t ctx) {
    fq_nmod_zero(result, ctx); 
    if (i == 0) return;

    fq_nmod_t gen;
    fq_nmod_init(gen, ctx);
    fq_nmod_gen(gen, ctx);

    fq_nmod_t power_of_a;
    fq_nmod_init(power_of_a, ctx);
    fq_nmod_one(power_of_a, ctx); 

    ulong p = fq_nmod_ctx_prime(ctx);
    ulong temp_i = i;

    while (temp_i > 0) {
        ulong remainder = temp_i % p;
        if (remainder != 0) {
            fq_nmod_t term;
            fq_nmod_init(term, ctx);
            fq_nmod_mul_ui(term, power_of_a, remainder, ctx);
            fq_nmod_add(result, result, term, ctx);
            fq_nmod_clear(term, ctx);
        }
        fq_nmod_mul(power_of_a, power_of_a, gen, ctx);
        temp_i /= p;
    }
    fq_nmod_clear(power_of_a, ctx);
    fq_nmod_clear(gen, ctx);
}

/**
 * @brief Função principal que orquestra a geração de todos os polinômios.
 */
fq_nmod_poly_t* generate_all_polynomials(long* poly_count, long max_degree, const fq_nmod_ctx_t ctx) {
    fmpz_t order_z, total_polys_z;
    fmpz_init(order_z);
    fmpz_init(total_polys_z);
    
    fq_nmod_ctx_order(order_z, ctx);
    long q = fmpz_get_si(order_z);

    // Calcula o número total de polinômios: q^(max_degree + 1)
    fmpz_pow_ui(total_polys_z, order_z, max_degree + 1);
    *poly_count = fmpz_get_si(total_polys_z);

    // Gera a lista de todos os elementos do corpo (para os coeficientes)
    fq_nmod_t* elements = (fq_nmod_t*) malloc(q * sizeof(fq_nmod_t));
    for (long i = 0; i < q; i++) {
        fq_nmod_init(elements[i], ctx);
        get_element_by_arithmetic(elements[i], i, ctx);
    }
    
    // Aloca memória para a lista de polinômios
    fq_nmod_poly_t* poly_list = (fq_nmod_poly_t*) malloc((*poly_count) * sizeof(fq_nmod_poly_t));
    for(long i=0; i < (*poly_count); i++) {
        fq_nmod_poly_init(poly_list[i], ctx);
    }

    // Prepara o pontapé inicial para a recursão
    fq_nmod_poly_t temp_poly;
    fq_nmod_poly_init(temp_poly, ctx);
    long start_index = 0;
    
    printf("Gerando %ld polinômios...\n", *poly_count);
    generate_recursive(poly_list, &start_index, temp_poly, 0, max_degree, elements, q, ctx);

    // Libera a memória temporária
    fq_nmod_poly_clear(temp_poly, ctx);
    for (long i = 0; i < q; i++) fq_nmod_clear(elements[i], ctx);
    free(elements);
    fmpz_clear(order_z);
    fmpz_clear(total_polys_z);
    
    return poly_list;
}

/**
 * @brief Função recursiva que constrói os polinômios coeficiente por coeficiente.
 */
void generate_recursive(
    fq_nmod_poly_t* poly_list, 
    long* current_index, 
    fq_nmod_poly_t current_poly, 
    long degree, 
    long max_degree, 
    const fq_nmod_t* elements, 
    long num_elements,
    const fq_nmod_ctx_t ctx) 
{
    // Caso base: se já definimos todos os coeficientes, o polinômio está pronto.
    if (degree > max_degree) {
        fq_nmod_poly_set(poly_list[*current_index], current_poly, ctx);
        (*current_index)++;
        return;
    }

    // Passo recursivo: para o grau atual, tente todos os elementos como coeficiente.
    for (long i = 0; i < num_elements; i++) {
        fq_nmod_poly_set_coeff(current_poly, degree, elements[i], ctx);
        generate_recursive(poly_list, current_index, current_poly, degree + 1, max_degree, elements, num_elements, ctx);
    }
}

/**
 * @brief Cria uma matriz 2D com os resultados da avaliação de polinômios.
 * @return Uma matriz alocada dinamicamente (ponteiro para ponteiro). O chamador DEVE liberar essa memória.
 */
fq_nmod_t** create_evaluation_matrix(long num_points, long num_polys, const fq_nmod_t* points, const fq_nmod_poly_t* polys, const fq_nmod_ctx_t ctx) {
    // Passo 1: Alocar o array de ponteiros (as linhas)
    fq_nmod_t** matrix = (fq_nmod_t**) malloc(num_points * sizeof(fq_nmod_t*));
    if (matrix == NULL) return NULL;

    // Passo 2: Para cada linha, alocar o array de elementos (as colunas) e inicializá-los
    for (long i = 0; i < num_points; i++) {
        matrix[i] = (fq_nmod_t*) malloc(num_polys * sizeof(fq_nmod_t));
        if (matrix[i] == NULL) { /* Lidar com erro de alocação */ return NULL; }

        for (long j = 0; j < num_polys; j++) {
            fq_nmod_init(matrix[i][j], ctx);
        }
    }

    // Passo 3: Preencher a matriz com as avaliações
    for (long i = 0; i < num_points; i++) {
        for (long j = 0; j < num_polys; j++) {
            fq_nmod_poly_evaluate_fq_nmod(matrix[i][j], polys[j], points[i], ctx);
        }
    }
    
    return matrix;
}

/**
 * @brief Libera a memória de uma matriz de avaliação criada com create_evaluation_matrix.
 */
void free_evaluation_matrix(fq_nmod_t** matrix, long num_points, long num_polys, const fq_nmod_ctx_t ctx) {
    if (matrix == NULL) return;

    for (long i = 0; i < num_points; i++) {
        if (matrix[i] != NULL) {
            // Primeiro, limpa cada elemento da linha
            for (long j = 0; j < num_polys; j++) {
                fq_nmod_clear(matrix[i][j], ctx);
            }
            // Depois, libera a própria linha
            free(matrix[i]);
        }
    }
    // Finalmente, libera o array de ponteiros
    free(matrix);
}