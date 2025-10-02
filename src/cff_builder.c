#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_poly.h"
#include "flint/nmod_poly.h"
#include "cff_builder.h"
#include "cff_file_generator.h"

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
void add_pair_to_list(element_pair** list, long* count, long* capacity, const fq_nmod_t x, const fq_nmod_t y, const fq_nmod_ctx_t ctx);
combination_partitions generate_combinations(const subfield_partition* partitions, int num_partitions, const fq_nmod_ctx_t ctx);
int** generate_single_cff(long* num_rows, const element_pair* combos,  long num_combos, const fq_nmod_t* points_for_eval,  long num_points, fq_nmod_t** evals, long num_polys, const fq_nmod_ctx_t ctx);
long find_element_index(const fq_nmod_t element, const fq_nmod_t* list, long list_count, const fq_nmod_ctx_t ctx);
void free_matrix(int** matrix, long rows);
generated_cffs generate_new_cff_blocks(long* Fq_steps, long* k_steps, int num_steps);
static void free_generated_cffs(generated_cffs* cffs);
void free_subfield_partitions(subfield_partition* partitions, int num_steps, const fq_nmod_ctx_t ctx);
void free_combination_partitions(combination_partitions* combos, const fq_nmod_ctx_t ctx);

void embeed_cff(long* Fq_steps, long* k_steps, int num_steps){
    int d0 = (Fq_steps[num_steps-2]-1)/(k_steps[num_steps-2]);
    long t0 = Fq_steps[num_steps-2] * Fq_steps[num_steps-2];
    long n0 = (long)pow(Fq_steps[num_steps-2], k_steps[num_steps-2] + 1);
    char filename0[100]; // buffer para a string
    snprintf(filename0, sizeof(filename0), "CFFs/%ld-CFF(%ld,%ld).txt", d0, t0, n0);

    long old_rows = 0, old_cols = 0;
    int** cff_old_old = read_cff_from_file(filename0, &old_rows, &old_cols);

    generated_cffs new_blocks = generate_new_cff_blocks(Fq_steps, k_steps, num_steps);

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

    int d1 = (Fq_steps[num_steps-1]-1)/(k_steps[num_steps-1]);
    long t1 = Fq_steps[num_steps-1] * Fq_steps[num_steps-1];
    long n1 = (long)pow(Fq_steps[num_steps-1], k_steps[num_steps-1] + 1);
    char filename1[100]; // buffer para a string
    snprintf(filename1, sizeof(filename1), "CFFs/%ld-CFF(%ld,%ld).txt", d1, t1, n1);

    write_cff_to_file(filename1, final_cff, new_total_rows, new_total_cols);

    // 5. LIMPEZA TOTAL DA MEMÓRIA
    printf("\nLimpando toda a memória.\n");
    free_matrix(cff_old_old, old_rows);
    free_generated_cffs(&new_blocks);
    free_matrix(final_cff, new_total_rows);
}


generated_cffs generate_new_cff_blocks(long* Fq_steps, long* k_steps, int num_steps) {
    generated_cffs result = {0};
    
    // --- ETAPA 0: PARÂMETROS E INICIALIZAÇÃO ---
    fmpz_t pz; fmpz_init(pz); 
    fmpz_set_ui(pz, Fq_steps[0]);

    nmod_poly_t mod_poly; 
    nmod_poly_init(mod_poly, Fq_steps[0]);
    nmod_poly_init(mod_poly, fmpz_get_ui(pz));
    nmod_poly_set_coeff_ui(mod_poly, 4, 1);
    nmod_poly_set_coeff_ui(mod_poly, 1, 1);
    nmod_poly_set_coeff_ui(mod_poly, 0, 1);

    fq_nmod_ctx_t ctx; 
    fq_nmod_ctx_init_modulus(ctx, mod_poly, "a");
    nmod_poly_clear(mod_poly);
    
    // --- ETAPA 1: PARTICIONAR ELEMENTOS ---
    subfield_partition* partitions = partition_by_subfields(Fq_steps, num_steps, ctx);

    // --- ETAPA 2: PARTICIONAR POLINÔMIOS ---
    fq_nmod_poly_t* old_polys = NULL;
    long num_old_polys_total = 0;

    for (int i = 0; i < num_steps - 1; i++) {
        long num_polys_in_step = 0;
        subfield_partition old_partition = partitions[i];

        fq_nmod_poly_t* old_polys_bySteps = generate_polynomials_from_coeffs(
            &num_polys_in_step, k_steps[i], old_partition.all_elements,
            old_partition.count_all, ctx);

        if (num_polys_in_step > 0) {
            fq_nmod_poly_t* unique_polys_in_step = (fq_nmod_poly_t*) malloc(num_polys_in_step * sizeof(fq_nmod_poly_t));
            if (unique_polys_in_step == NULL) { exit(EXIT_FAILURE); }
            long num_unique_in_step = 0;

            for (long k = 0; k < num_polys_in_step; k++) {
                fq_nmod_poly_init(unique_polys_in_step[k], ctx);
            }

            for (long j = 0; j < num_polys_in_step; j++) {
                if (!fq_nmod_poly_is_in_list(old_polys_bySteps[j], old_polys, num_old_polys_total, ctx)) {
                    
                    fq_nmod_poly_set(unique_polys_in_step[num_unique_in_step], old_polys_bySteps[j], ctx);
                    num_unique_in_step++;
                }
            }

            if (num_unique_in_step > 0) {
                fq_nmod_poly_t* temp = realloc(old_polys, (num_old_polys_total + num_unique_in_step) * sizeof(fq_nmod_poly_t));
                if (temp == NULL) {
                    fprintf(stderr, "Erro ao realocar memória!\n");
                    for (long k = 0; k < num_polys_in_step; k++) fq_nmod_poly_clear(unique_polys_in_step[k], ctx);
                    free(unique_polys_in_step);
                    exit(EXIT_FAILURE);
                }
                old_polys = temp;

                for (long j = 0; j < num_unique_in_step; j++) {
                    long dest_index = num_old_polys_total + j;
                    fq_nmod_poly_init(old_polys[dest_index], ctx);
                    fq_nmod_poly_set(old_polys[dest_index], unique_polys_in_step[j], ctx);
                }

                num_old_polys_total += num_unique_in_step;
            }

            for (long k = 0; k < num_polys_in_step; k++) {
                fq_nmod_poly_clear(unique_polys_in_step[k], ctx);
            }
            free(unique_polys_in_step);
        }

        for (long j = 0; j < num_polys_in_step; j++) {
            fq_nmod_poly_clear(old_polys_bySteps[j], ctx);
        }
        free(old_polys_bySteps);
    }


    subfield_partition all_partition = partitions[num_steps - 1];
    long num_all_polys = 0;
    fq_nmod_poly_t* all_polys = generate_polynomials_from_coeffs(&num_all_polys, k_steps[num_steps-1], all_partition.all_elements, all_partition.count_all, ctx);
    
    fq_nmod_poly_t* new_polys = (fq_nmod_poly_t*) malloc(num_all_polys * sizeof(fq_nmod_poly_t));
    long num_new_polys = 0;
    for (long i = 0; i < num_all_polys; i++) {
        if (!fq_nmod_poly_is_in_list(all_polys[i], old_polys, num_old_polys_total, ctx)) {
            fq_nmod_poly_init(new_polys[num_new_polys], ctx);
            fq_nmod_poly_set(new_polys[num_new_polys], all_polys[i], ctx);
            num_new_polys++;
        }
    }

    /*
    for(long i = 0; i < num_old_polys_total; i++){
        fq_nmod_poly_print_pretty(old_polys[i], "x", ctx);
        printf("\n");
    }

    for(long i = 0; i < num_new_polys; i++){
        fq_nmod_poly_print_pretty(new_polys[i], "x", ctx);
        printf("\n");
    }

    printf("%ld\n", num_new_polys);
    */

    // --- ETAPA 3: GERAR COMBINAÇÕES DE ELEMENTOS ---
    combination_partitions combos = generate_combinations(partitions, num_steps, ctx);

    /*
    printf("--- Combos Old (Total: %ld) ---\n", combos.count_old);
    for (long i = 0; i < combos.count_old; i++) {
        printf("Par %ld: (", i);
        fq_nmod_print_pretty(combos.combos_old[i].x, ctx); // Correto: .x
        printf(", ");
        fq_nmod_print_pretty(combos.combos_old[i].y, ctx); // Correto: .y
        printf(")\n");
    }
    printf("\n");

    printf("--- Combos New (Total: %ld) ---\n", combos.count_new);
    for (long i = 0; i < combos.count_new; i++) {
        printf("Par %ld: (", i);
        fq_nmod_print_pretty(combos.combos_new[i].x, ctx); // Correto: combos_new e .x
        printf(", ");
        fq_nmod_print_pretty(combos.combos_new[i].y, ctx); // Correto: combos_new e .y
        printf(")\n");
    }
    printf("----------------------------------\n\n");
    */

    // --- ETAPA 4: CRIAR MATRIZES DE AVALIAÇÃO ---
    fq_nmod_t* points_for_eval = all_partition.all_elements;
    long num_points = all_partition.count_all;
    fq_nmod_t** evals_old = create_evaluation_matrix(num_points, num_old_polys_total, points_for_eval, old_polys, ctx);
    fq_nmod_t** evals_new = create_evaluation_matrix(num_points, num_new_polys, points_for_eval, new_polys, ctx);

    // --- ETAPA 5: GERAR MATRIZES CFF FINAIS ---
    result.cff_old_new = generate_single_cff(&result.rows_old_new, combos.combos_old, combos.count_old, points_for_eval, num_points, evals_new, num_new_polys, ctx);
    result.cols_old_new = num_new_polys;
    
    result.cff_new_old = generate_single_cff(&result.rows_new_old, combos.combos_new, combos.count_new, points_for_eval, num_points, evals_old, num_old_polys_total, ctx);
    result.cols_new_old = num_old_polys_total;

    result.cff_new = generate_single_cff(&result.rows_new, combos.combos_new, combos.count_new, points_for_eval, num_points, evals_new, num_new_polys, ctx);
    result.cols_new = num_new_polys;
    
    // --- ETAPA 6: LIMPEZA ---
    free_evaluation_matrix(evals_old, num_points, num_old_polys_total, ctx);
    free_evaluation_matrix(evals_new, num_points, num_new_polys, ctx);
    free_combination_partitions(&combos, ctx);
    if (old_polys != NULL) { for (long i = 0; i < num_old_polys_total; i++) fq_nmod_poly_clear(old_polys[i], ctx); free(old_polys); }
    if (new_polys != NULL) { for (long i = 0; i < num_new_polys; i++) fq_nmod_poly_clear(new_polys[i], ctx); free(new_polys); }
    if (all_polys != NULL) { for (long i = 0; i < num_all_polys; i++) fq_nmod_poly_clear(all_polys[i], ctx); free(all_polys); }
    free_subfield_partitions(partitions, num_steps, ctx);
    fmpz_clear(pz);
    fq_nmod_ctx_clear(ctx);

    return result;
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
 * 
 *  ----- REVER REALLOC -----
 * 
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

/**
 * @brief Gera listas de pares 'old' e 'new' com base nas partições de subcorpos.
 */
combination_partitions generate_combinations(const subfield_partition* partitions, int num_partitions, const fq_nmod_ctx_t ctx) {
    combination_partitions result = {0}; // Inicializa tudo com NULL/0

    fq_nmod_t* all_accumulated_elements = NULL;
    long accumulated_count = 0;
    long accumulated_capacity = 0;

    for (int i = 0; i < num_partitions; i++) {
        fq_nmod_t* current_only_elements = partitions[i].only_elements;
        long current_only_count = partitions[i].count_only;

        // Decide se os novos pares vão para a lista 'old' ou 'new'
        element_pair** target_list = (i == num_partitions - 1) ? &result.combos_new : &result.combos_old;
        long* target_count = (i == num_partitions - 1) ? &result.count_new : &result.count_old;
        long* target_capacity = (i == num_partitions - 1) ? &result.capacity_new : &result.capacity_old;

        if (i == 0) {
            // Caso base: pares dentro do menor subcorpo
            for (long ix = 0; ix < current_only_count; ix++) {
                for (long iy = 0; iy < current_only_count; iy++) {
                    add_pair_to_list(target_list, target_count, target_capacity, current_only_elements[ix], current_only_elements[iy], ctx);
                }
            }
        } else {
            // Casos intermediário e final: pares entre o acumulado e o novo
            // Pares (acumulado, novo)
            for (long ix = 0; ix < accumulated_count; ix++) {
                for (long iy = 0; iy < current_only_count; iy++) {
                    add_pair_to_list(target_list, target_count, target_capacity, all_accumulated_elements[ix], current_only_elements[iy], ctx);
                }
            }
            // Pares (novo, acumulado) - ATENÇÃO: o acumulado precisa ser o já estendido
            // Primeiro, estendemos a lista acumulada
            for (long j = 0; j < current_only_count; j++) {
                add_element_to_list(&all_accumulated_elements, &accumulated_count, &accumulated_capacity, current_only_elements[j], ctx);
            }
            // Agora, fazemos os pares com a lista estendida
            for (long ix = 0; ix < current_only_count; ix++) {
                for (long iy = 0; iy < accumulated_count; iy++) {
                     add_pair_to_list(target_list, target_count, target_capacity, current_only_elements[ix], all_accumulated_elements[iy], ctx);
                }
            }
        }
        
        // Se for o caso base, a lista acumulada ainda não foi preenchida
        if (i == 0) {
             for (long j = 0; j < current_only_count; j++) {
                add_element_to_list(&all_accumulated_elements, &accumulated_count, &accumulated_capacity, current_only_elements[j], ctx);
            }
        }
    }
    
    // Libera a lista acumulada temporária
    for (long i = 0; i < accumulated_count; i++) {
        fq_nmod_clear(all_accumulated_elements[i], ctx);
    }
    free(all_accumulated_elements);

    return result;
}

/**
 * @brief Adiciona um par de elementos (x, y) a uma lista dinâmica de pares.
 */
void add_pair_to_list(element_pair** list, long* count, long* capacity, const fq_nmod_t x, const fq_nmod_t y, const fq_nmod_ctx_t ctx) {
    if (*count >= *capacity) {
        *capacity = (*capacity == 0) ? 8 : (*capacity) * 2;
        *list = realloc(*list, (*capacity) * sizeof(element_pair));
        for (long i = *count; i < *capacity; i++) {
            fq_nmod_init((*list)[i].x, ctx);
            fq_nmod_init((*list)[i].y, ctx);
        }
    }
    fq_nmod_set((*list)[*count].x, x, ctx);
    fq_nmod_set((*list)[*count].y, y, ctx);
    (*count)++;
}

/**
 * @brief Gera uma matriz CFF (binária) a partir de combinações e avaliações.
 */
int** generate_single_cff(long* num_rows, const element_pair* combos,  long num_combos, const fq_nmod_t* points_for_eval,  long num_points, fq_nmod_t** evals, long num_polys, const fq_nmod_ctx_t ctx) {
    *num_rows = num_combos;
    if (num_combos == 0) return NULL;

    // Aloca a matriz CFF binária
    int** cff_matrix = (int**) malloc(num_combos * sizeof(int*));

    fq_nmod_t x, y;
    fq_nmod_init(x, ctx);
    fq_nmod_init(y, ctx);

    // Itera por cada par (x,y)
    for (long i = 0; i < num_combos; i++) {
        fq_nmod_set(x, combos[i].x, ctx);
        fq_nmod_set(y, combos[i].y, ctx);

        // Encontra o índice da linha na matriz de avaliação que corresponde a 'x'
        long row_idx = find_element_index(x, points_for_eval, num_points, ctx);
        
        // Aloca a linha na nossa matriz CFF de saída
        cff_matrix[i] = (int*) malloc(num_polys * sizeof(int));

        if (row_idx != -1) {
            // Se encontramos o elemento, comparamos a linha de avaliação com 'y'
            fq_nmod_t* eval_row = evals[row_idx];
            for (long j = 0; j < num_polys; j++) {
                if (fq_nmod_equal(eval_row[j], y, ctx)) {
                    cff_matrix[i][j] = 1;
                } else {
                    cff_matrix[i][j] = 0;
                }
            }
        } else {
            // Se 'x' não estava na lista de pontos avaliados (não deveria acontecer),
            // a linha será toda de zeros.
            for (long j = 0; j < num_polys; j++) {
                cff_matrix[i][j] = 0;
            }
        }
    }

    fq_nmod_clear(x, ctx);
    fq_nmod_clear(y, ctx);

    return cff_matrix;
}

/**
 * @brief Encontra o índice de um elemento em uma lista de elementos.
 * @return O índice se encontrado, ou -1 se não encontrado.
 */
long find_element_index(const fq_nmod_t element, const fq_nmod_t* list, long list_count, const fq_nmod_ctx_t ctx) {
    for (long i = 0; i < list_count; i++) {
        if (fq_nmod_equal(element, list[i], ctx)) {
            return i; // Encontrado!
        }
    }
    return -1; // Não encontrado
}

void free_subfield_partitions(subfield_partition* partitions, int num_steps, const fq_nmod_ctx_t ctx) {
    if (!partitions) return;
    for (int i = 0; i < num_steps; i++) {
        for (long j = 0; j < partitions[i].capacity_all; j++) fq_nmod_clear(partitions[i].all_elements[j], ctx);
        for (long j = 0; j < partitions[i].capacity_only; j++) fq_nmod_clear(partitions[i].only_elements[j], ctx);
        free(partitions[i].all_elements);
        free(partitions[i].only_elements);
    }
    free(partitions);
}

void free_combination_partitions(combination_partitions* combos, const fq_nmod_ctx_t ctx) {
    if (!combos) return;
    for (long i = 0; i < combos->capacity_old; i++) { fq_nmod_clear(combos->combos_old[i].x, ctx); fq_nmod_clear(combos->combos_old[i].y, ctx); }
    free(combos->combos_old);
    for (long i = 0; i < combos->capacity_new; i++) { fq_nmod_clear(combos->combos_new[i].x, ctx); fq_nmod_clear(combos->combos_new[i].y, ctx); }
    free(combos->combos_new);
}

void free_generated_cffs(generated_cffs* cffs) {
    if (!cffs) return;
    if (cffs->cff_old_new) { for (long i = 0; i < cffs->rows_old_new; i++) free(cffs->cff_old_new[i]); free(cffs->cff_old_new); }
    if (cffs->cff_new_old) { for (long i = 0; i < cffs->rows_new_old; i++) free(cffs->cff_new_old[i]); free(cffs->cff_new_old); }
    if (cffs->cff_new) { for (long i = 0; i < cffs->rows_new; i++) free(cffs->cff_new[i]); free(cffs->cff_new); }
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

void free_matrix(int** matrix, long rows) {
    if (!matrix) return;
    for (long i = 0; i < rows; i++) free(matrix[i]);
    free(matrix);
}