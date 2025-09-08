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
void add_poly_to_list(fq_nmod_poly_t** list, long* count, long* capacity, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx);
subfield_partition* partition_by_subfields(const long* Fq_steps, int num_steps, const fq_nmod_ctx_t ctx);
fq_nmod_poly_t* generate_polynomials_from_coeffs(long* poly_count, long max_degree, const fq_nmod_t* coeffs, long num_coeffs, const fq_nmod_ctx_t ctx);
void generate_recursive_sorted(fq_nmod_poly_t* poly_list, long* current_index, fq_nmod_poly_t current_poly, long degree, const fq_nmod_t* elements, long num_elements, const fq_nmod_ctx_t ctx);
fq_nmod_t** create_evaluation_matrix(long num_points, long num_polys, const fq_nmod_t* points, const fq_nmod_poly_t* polys, const fq_nmod_ctx_t ctx);
void free_evaluation_matrix(fq_nmod_t** matrix, long num_points, long num_polys, const fq_nmod_ctx_t ctx);
int fq_nmod_poly_is_in_list(const fq_nmod_poly_t poly, const fq_nmod_poly_t* list, long list_count, const fq_nmod_ctx_t ctx);


int main() {
    // --- Parâmetros ---
    fmpz_t pz; fmpz_init(pz); fmpz_set_ui(pz, 2);
    long D = 2; // F_q^d
    long Fq_steps[] = {2,4}; 
    int num_steps = 2;

    // --- Inicialização do Corpo ---
    nmod_poly_t mod_poly; nmod_poly_init(mod_poly, 2);
    nmod_poly_set_coeff_ui(mod_poly, 2, 1); 
    nmod_poly_set_coeff_ui(mod_poly, 1, 1); 
    nmod_poly_set_coeff_ui(mod_poly, 0, 1);
    fq_nmod_ctx_t ctx; fq_nmod_ctx_init_modulus(ctx, mod_poly, "a");
    printf("Contexto inicializado para F_%ld^%ld.\n\n", fmpz_get_ui(pz), D);
    nmod_poly_clear(mod_poly);
    
    // --- Particionamento de Elementos ---
    subfield_partition* partitions = partition_by_subfields(Fq_steps, num_steps, ctx);

    // --- GERAÇÃO DE POLINÔMIOS ---
    long max_poly_degree = 1;
    
    // Inicializa uma lista dinâmica para 'all_polys'
    fq_nmod_poly_t* all_polys = NULL;
    long all_polys_count = 0;
    long all_polys_capacity = 0;

    printf("--- Gerando Polinômios Particionados ---\n");
    // Itera sobre as partições de subcorpos (F_q, F_qˆd, ...)
    for(int i = 0; i < num_steps; i++) {
        long num_polys_level = 0;
        subfield_partition current_partition = partitions[i];
        
        // Gera todos os polinômios com coeficientes do subcorpo atual
        fq_nmod_poly_t* level_polys = generate_polynomials_from_coeffs(
            &num_polys_level, max_poly_degree, 
            current_partition.all_elements, current_partition.count_all, ctx
        );

        // Adiciona os polinômios gerados à lista principal, se ainda não estiverem lá
        for(long j = 0; j < num_polys_level; j++) {
            int found = fq_nmod_poly_is_in_list(level_polys[j], all_polys, all_polys_count, ctx);
            if (!found) {
                add_poly_to_list(&all_polys, &all_polys_count, &all_polys_capacity, level_polys[j], ctx);
            } 
        }

        // Libera a lista temporária de polinômios deste nível
        for (long j = 0; j < num_polys_level; j++) {
            fq_nmod_poly_clear(level_polys[j], ctx);
        }
        free(level_polys);
    }

    // --- Impressão dos Polinômios Únicos Gerados ---
    if (all_polys != NULL) {
        printf("\n--- %ld Polinômios Únicos Gerados (Grau <= %ld) ---\n", all_polys_count, max_poly_degree);
        for (long i = 0; i < all_polys_count; i++) {
            printf("  P_%ld(x) = ", i);
            fq_nmod_poly_print_pretty(all_polys[i], "x", ctx);
            printf("\n");
        }
    }

    // --- AVALIAÇÃO DE POLINÔMIOS ---
    if (all_polys != NULL && partitions != NULL && num_steps > 1) {
        printf("\n--- Avaliando Polinômios ---\n");
        subfield_partition f4_partition = partitions[num_steps - 1]; 
        long num_points_to_eval = f4_partition.count_all;

        fq_nmod_t** evals = create_evaluation_matrix(
            num_points_to_eval, all_polys_count, 
            f4_partition.all_elements, all_polys, ctx
        );

        if (evals != NULL) {
            printf("Matriz de Avaliação Completa [ponto][polinômio]:\n");
            
            for (long i = 0; i < num_points_to_eval; i++) {
                printf("  Ponto "); fq_nmod_print_pretty(f4_partition.all_elements[i], ctx); printf(": [ ");
                for (long j = 0; j < all_polys_count; j++) {
                    fq_nmod_print_pretty(evals[i][j], ctx);
                    if (j < all_polys_count - 1) printf(", ");
                }
                printf(" ]\n");
            }
            free_evaluation_matrix(evals, num_points_to_eval, all_polys_count, ctx);
        }
    }
    
    // --- Limpeza de Memória ---
    printf("\nLimpando a memória...\n");
    if (all_polys != NULL) {
        for (long i = 0; i < all_polys_count; i++) {
            fq_nmod_poly_clear(all_polys[i], ctx);
        }
        free(all_polys);
    }
    if (partitions != NULL) {
        for (int i = 0; i < num_steps; i++) {
            for (long j = 0; j < partitions[i].capacity_all; j++) fq_nmod_clear(partitions[i].all_elements[j], ctx);
            for (long j = 0; j < partitions[i].capacity_only; j++) fq_nmod_clear(partitions[i].only_elements[j], ctx);
            free(partitions[i].all_elements); free(partitions[i].only_elements);
        }
        free(partitions);
    }
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
fq_nmod_poly_t* generate_polynomials_from_coeffs(long* poly_count, long max_degree, const fq_nmod_t* coeffs, long num_coeffs, const fq_nmod_ctx_t ctx) {
    fmpz_t num_coeffs_z, total_polys_z;
    fmpz_init(num_coeffs_z);
    fmpz_init(total_polys_z);
    fmpz_set_si(num_coeffs_z, num_coeffs);

    // Calcula o número total de polinômios: (num_coeffs)^(max_degree + 1)
    fmpz_pow_ui(total_polys_z, num_coeffs_z, max_degree + 1);
    *poly_count = fmpz_get_si(total_polys_z);

    // Aloca memória para a lista de polinômios
    fq_nmod_poly_t* poly_list = (fq_nmod_poly_t*) malloc((*poly_count) * sizeof(fq_nmod_poly_t));
    for(long i=0; i < (*poly_count); i++) {
        fq_nmod_poly_init(poly_list[i], ctx);
    }

    fq_nmod_poly_t temp_poly;
    fq_nmod_poly_init(temp_poly, ctx);
    long start_index = 0;
    
    printf("Gerando %ld polinômios com %ld coeficientes possíveis...\n", *poly_count, num_coeffs);
    generate_recursive_sorted(poly_list, &start_index, temp_poly, max_degree, coeffs, num_coeffs, ctx);

    fq_nmod_poly_clear(temp_poly, ctx);
    fmpz_clear(num_coeffs_z);
    fmpz_clear(total_polys_z);
    
    return poly_list;
}

/**
 * @brief Função recursiva que constrói os polinômios, priorizando o termo de maior grau.
 * Isso resulta em uma lista de polinômios mais ordenada e intuitiva.
 */
void generate_recursive_sorted(
    fq_nmod_poly_t* poly_list, 
    long* current_index, 
    fq_nmod_poly_t current_poly, 
    long degree, // Começará com max_degree e diminuirá
    const fq_nmod_t* elements, 
    long num_elements,
    const fq_nmod_ctx_t ctx) 
{
    // Caso base: se já definimos todos os coeficientes (de max_degree até 0), o polinômio está pronto.
    if (degree < 0) {
        fq_nmod_poly_set(poly_list[*current_index], current_poly, ctx);
        (*current_index)++;
        return;
    }

    // Passo recursivo: para o grau atual, tente todos os elementos como coeficiente.
    for (long i = 0; i < num_elements; i++) {
        fq_nmod_poly_set_coeff(current_poly, degree, elements[i], ctx);
        // Chama a recursão para o próximo grau (um grau menor).
        generate_recursive_sorted(poly_list, current_index, current_poly, degree - 1, elements, num_elements, ctx);
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

/**
 * @brief Verifica se um polinômio existe em uma lista de polinômios.
 * @return 1 se encontrado, 0 caso contrário.
 */
int fq_nmod_poly_is_in_list(const fq_nmod_poly_t poly, const fq_nmod_poly_t* list, long list_count, const fq_nmod_ctx_t ctx) {
    for (long i = 0; i < list_count; i++) {
        if (fq_nmod_poly_equal(poly, list[i], ctx)) {
            return 1; // Encontrado
        }
    }
    return 0; // Não encontrado
}

/**
 * @brief Adiciona um polinômio a um array dinâmico de polinômios.
 */
void add_poly_to_list(fq_nmod_poly_t** list, long* count, long* capacity, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx) {
    if (*count >= *capacity) {
        *capacity = (*capacity == 0) ? 8 : (*capacity) * 2;
        *list = realloc(*list, (*capacity) * sizeof(fq_nmod_poly_t));
        for (long i = *count; i < *capacity; i++) {
            fq_nmod_poly_init((*list)[i], ctx);
        }
    }
    fq_nmod_poly_set((*list)[*count], poly, ctx);
    (*count)++;
}