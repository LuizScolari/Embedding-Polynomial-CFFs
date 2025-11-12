#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <glib.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_poly.h"
#include "flint/nmod_poly.h"
#include "cff_builder.h"
#include "cff_file_generator.h"

static fq_nmod_ctx_t global_ctx;
static int global_ctx_initialized = 0;

// --- Protótipos de Funções ---
void get_element_by_arithmetic(fq_nmod_t result, ulong i, fq_nmod_ctx_t ctx);
void add_element_to_list(fq_nmod_t** list, long* count, long* capacity, const fq_nmod_t element, const fq_nmod_ctx_t ctx);
void add_poly_to_list(fq_nmod_poly_t** list, long* count, long* capacity, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx);
subfield_partition* partition_by_subfields(const long* Fq_steps, int num_steps, const fq_nmod_ctx_t ctx);
fq_nmod_poly_t* generate_polynomials_from_coeffs(long* poly_count, long max_degree, const fq_nmod_t* coeffs, long num_coeffs, const fq_nmod_ctx_t ctx);
void generate_recursive_sorted(fq_nmod_poly_t* poly_list, long* current_index, fq_nmod_poly_t current_poly, long degree, const fq_nmod_t* elements, long num_elements, const fq_nmod_ctx_t ctx);
int fq_nmod_poly_is_in_list(const fq_nmod_poly_t poly, const fq_nmod_poly_t* list, long list_count, const fq_nmod_ctx_t ctx);
void add_pair_to_list(element_pair** list, long* count, long* capacity, const fq_nmod_t x, const fq_nmod_t y, const fq_nmod_ctx_t ctx);
combination_partitions generate_combinations(const subfield_partition* partitions, int num_partitions, const fq_nmod_ctx_t ctx);
int** generate_single_cff(long* num_rows, const element_pair* combos,  long num_combos, GHashTable* inverted_evals, long num_polys, const fq_nmod_ctx_t ctx);
long find_element_index(const fq_nmod_t element, const fq_nmod_t* list, long list_count, const fq_nmod_ctx_t ctx);
void free_matrix(int** matrix, long rows);
generated_cffs generate_new_cff_blocks(long* Fq_steps, long* k_steps, int num_steps);
static void free_generated_cffs(generated_cffs* cffs);
void free_subfield_partitions(subfield_partition* partitions, int num_steps, const fq_nmod_ctx_t ctx);
void free_combination_partitions(combination_partitions* combos, const fq_nmod_ctx_t ctx);
GHashTable* create_inverted_evaluation_index(long num_points, long num_polys, const fq_nmod_t* points, const fq_nmod_poly_t* polys, const fq_nmod_ctx_t ctx);
guint fq_nmod_hash_func(gconstpointer key);
gboolean fq_nmod_equal_func(gconstpointer a, gconstpointer b);
void g_array_destroy_wrapper(gpointer data);
void g_hash_table_destroy_wrapper(gpointer data);

void generate_cff(char construction, long fq, long k) {
    int d0 = (fq - 1) / k;
    long t0 = fq * fq;
    long n0 = (long)pow(fq, k + 1);
    char filename0[100];
    snprintf(filename0, sizeof(filename0), "CFFs/%d-CFF(%ld,%ld).txt", d0, t0, n0);
    printf("Gerando CFF inicial em '%s'...\n", filename0);

    long fq_array[1] = { fq };
    long k_array[1] = { k };
    int num_steps = 1;

    generated_cffs new_blocks = generate_new_cff_blocks(fq_array, k_array, num_steps);

    int** final_cff = new_blocks.cff_new;
    long final_rows = new_blocks.rows_new;
    long final_cols = new_blocks.cols_new;

    if (final_cff == NULL) {
        printf("Erro: Falha ao gerar a matriz CFF inicial.\n");
        return;
    }
    printf("Matriz CFF inicial de %ldx%ld gerada.\n", final_rows, final_cols);

    write_cff_to_file(filename0, construction, fq_array, 1, k_array, 1, final_cff, final_rows, final_cols);

    new_blocks.cff_old_new = NULL;
    new_blocks.cff_new_old = NULL;

    free_generated_cffs(&new_blocks); 
}

void embeed_cff(char construction, long* Fq_steps, long* k_steps){
    int d0 = (Fq_steps[0]-1)/(k_steps[0]);
    long t0 = Fq_steps[0] * Fq_steps[0];
    long n0 = (long)pow(Fq_steps[0], k_steps[0] + 1);
    char filename0[100]; // buffer para a string
    snprintf(filename0, sizeof(filename0), "CFFs/%d-CFF(%ld,%ld).txt", d0, t0, n0);

    long old_rows = 0, old_cols = 0;
    struct cff_parameters* params = read_parameters(filename0);
    if (params == NULL) {
        printf("Erro ao ler parâmetros do arquivo %s\n", filename0);
        return;
    }

    int** cff_old_old = read_cff_from_file(filename0, &old_rows, &old_cols);

    int new_fqs_count = params->fqs_count + 1;
    int new_ks_count = params->ks_count + 1;
    long* new_Fq_steps = (long*) malloc(new_fqs_count * sizeof(long));
    long* new_k_steps  = (long*) malloc(new_ks_count * sizeof(long));

    if (new_Fq_steps == NULL || new_k_steps == NULL) {
        printf("Erro: Falha ao alocar memória para conversão dos passos.\n");
        free(params->Fqs);
        free(params->ks);
        free(params);
        free(new_Fq_steps); 
        free(new_k_steps);
        free_matrix(cff_old_old, old_rows);
        return;
    }

    for (int i = 0; i < new_fqs_count; i++) {
        if (i < params->fqs_count) { 
            new_Fq_steps[i] = (long)params->Fqs[i];
        } else {
            new_Fq_steps[i] = Fq_steps[1];
        }
    }


    for (int i = 0; i < new_ks_count; i++) {
        if (i < params->ks_count) { 
            new_k_steps[i] = (long)params->ks[i];
        } else {
            new_k_steps[i] = k_steps[1];
        }
    }

    generated_cffs new_blocks = generate_new_cff_blocks(new_Fq_steps, new_k_steps, new_fqs_count);

    long new_total_rows = old_rows + new_blocks.rows_new_old;

    long width_top = old_cols + new_blocks.cols_old_new;
    long width_bottom = new_blocks.cols_new_old + new_blocks.cols_new;
    long new_total_cols = (width_top > width_bottom) ? width_top : width_bottom;
    
    int** final_cff = (int**) malloc(new_total_rows * sizeof(int*));
    for (long i = 0; i < new_total_rows; i++) {
        final_cff[i] = (int*) calloc(new_total_cols, sizeof(int));
    }

    // Copia bloco 1: CFF_old_old
    for(long i=0; i < old_rows; i++) memcpy(final_cff[i], cff_old_old[i], old_cols * sizeof(int));

    // Copia bloco 2: CFF_old_new
    for(long i=0; i < new_blocks.rows_old_new; i++) memcpy(&final_cff[i][old_cols], new_blocks.cff_old_new[i], new_blocks.cols_old_new * sizeof(int));
    
    // Copia bloco 3: CFF_new_old
    for(long i=0; i < new_blocks.rows_new_old; i++) memcpy(final_cff[old_rows + i], new_blocks.cff_new_old[i], new_blocks.cols_new_old * sizeof(int));
    
    // Copia bloco 4: CFF_new
    for(long i=0; i < new_blocks.rows_new; i++) memcpy(&final_cff[old_rows + i][new_blocks.cols_new_old], new_blocks.cff_new[i], new_blocks.cols_new * sizeof(int));

    int d1 = (Fq_steps[1]-1)/(k_steps[1]);
    long t1 = Fq_steps[1] * Fq_steps[1];
    long n1 = (long)pow(Fq_steps[1], k_steps[1] + 1);
    char filename1[100]; // buffer para a string
    snprintf(filename1, sizeof(filename1), "CFFs/%d-CFF(%ld,%ld).txt", d1, t1, n1);

    write_cff_to_file(filename1, params->construction, new_Fq_steps, new_fqs_count, new_k_steps, new_ks_count, final_cff, new_total_rows, new_total_cols);

    // 5. LIMPEZA TOTAL DA MEMÓRIA
    free_matrix(cff_old_old, old_rows);
    free_generated_cffs(&new_blocks);
    free_matrix(final_cff, new_total_rows);
    free(new_Fq_steps);
    free(new_k_steps);
    free(params->Fqs);
    free(params->ks);
    free(params);
}

generated_cffs generate_new_cff_blocks(long* Fq_steps, long* k_steps, int num_steps) {
    generated_cffs result = {0};
    
    // --- ETAPA 0: PARÂMETROS E INICIALIZAÇÃO ---
    ulong p = Fq_steps[0];
    long q = Fq_steps[num_steps-1]; 
    long n = (long)round(log(q) / log(p));

    fmpz_t pz;
    fmpz_init(pz);
    fmpz_set_ui(pz, p);

    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx, p, n, "a");
    
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
    GHashTable* inverted_index_old = create_inverted_evaluation_index(num_points, num_old_polys_total, points_for_eval, old_polys, ctx);
    GHashTable* inverted_index_new = create_inverted_evaluation_index(num_points, num_new_polys, points_for_eval, new_polys, ctx);

    // O contexto global para as funções de comparação da GHashTable ainda é necessário.
    if (!global_ctx_initialized) {
        fq_nmod_ctx_init_modulus(global_ctx, fq_nmod_ctx_modulus(ctx), "a");
        global_ctx_initialized = 1;
    }

    // --- ETAPA 5: GERAR MATRIZES CFF FINAIS ---
    result.cff_old_new = generate_single_cff(&result.rows_old_new, combos.combos_old, combos.count_old, inverted_index_new, num_new_polys, ctx);
    result.cols_old_new = num_new_polys;
    
    result.cff_new_old = generate_single_cff(&result.rows_new_old, combos.combos_new, combos.count_new, inverted_index_old, num_old_polys_total, ctx);
    result.cols_new_old = num_old_polys_total;

    result.cff_new = generate_single_cff(&result.rows_new, combos.combos_new, combos.count_new, inverted_index_new, num_new_polys, ctx);
    result.cols_new = num_new_polys;
    
    // --- ETAPA 6: LIMPEZA ---
    g_hash_table_destroy(inverted_index_old);
    g_hash_table_destroy(inverted_index_new);
    free_combination_partitions(&combos, ctx);
    if (old_polys != NULL) { for (long i = 0; i < num_old_polys_total; i++) fq_nmod_poly_clear(old_polys[i], ctx); free(old_polys); }
    if (new_polys != NULL) { for (long i = 0; i < num_new_polys; i++) fq_nmod_poly_clear(new_polys[i], ctx); free(new_polys); }
    if (all_polys != NULL) { for (long i = 0; i < num_all_polys; i++) fq_nmod_poly_clear(all_polys[i], ctx); free(all_polys); }
    free_subfield_partitions(partitions, num_steps, ctx);
    fq_nmod_ctx_clear(ctx);
    fq_nmod_ctx_clear(global_ctx);
    global_ctx_initialized = 0;

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
int** generate_single_cff(long* num_rows, const element_pair* combos,  long num_combos, GHashTable* inverted_evals, long num_polys, const fq_nmod_ctx_t ctx) {
    *num_rows = num_combos;
    if (num_combos == 0) return NULL;

    int** cff_matrix = (int**) malloc(num_combos * sizeof(int*));
    if (cff_matrix == NULL) exit(EXIT_FAILURE);

    // Itera por cada par (x,y)
    for (long i = 0; i < num_combos; i++) {
        const fq_nmod_t* x = &combos[i].x;
        const fq_nmod_t* y = &combos[i].y;

        // Aloca a linha na nossa matriz CFF e inicializa com zeros.
        // calloc é perfeito para isso.
        cff_matrix[i] = (int*) calloc(num_polys, sizeof(int));
        if (cff_matrix[i] == NULL) exit(EXIT_FAILURE);

        // O(1) na média para encontrar a tabela interna para 'x'
        GHashTable* inner_hash = g_hash_table_lookup(inverted_evals, x);

        if (inner_hash != NULL) {
            // O(1) na média para encontrar a lista de índices para 'y'
            GArray* indices = g_hash_table_lookup(inner_hash, y);

            if (indices != NULL) {
                // Agora, em vez de percorrer todos os 'num_polys',
                // percorremos apenas a pequena lista de índices corretos.
                for (guint k = 0; k < indices->len; k++) {
                    long j = g_array_index(indices, long, k);
                    cff_matrix[i][j] = 1;
                }
            }
        }
        // Se inner_hash ou indices for NULL, a linha já está corretamente preenchida com zeros
        // graças ao calloc, então não precisamos de um 'else'.
    }

    return cff_matrix;
}

// Função para criar o índice invertido
GHashTable* create_inverted_evaluation_index(long num_points, long num_polys, const fq_nmod_t* points, const fq_nmod_poly_t* polys, const fq_nmod_ctx_t ctx) {
    
    // Tabela externa: mapeia x -> (tabela interna)
    // Precisamos de destrutores personalizados para as chaves e valores
    GHashTable* inverted_index = g_hash_table_new_full(
        fq_nmod_hash_func, 
        fq_nmod_equal_func, 
        g_free, // Destruidor para as chaves (que serão cópias)
        g_hash_table_destroy_wrapper // Destruidor para os valores (as tabelas internas)
    );

    fq_nmod_t y_eval;
    fq_nmod_init(y_eval, ctx);

    for (long i = 0; i < num_points; i++) {
        const fq_nmod_t* x = &points[i];

        // Tabela interna: mapeia y -> (lista de índices j)
        GHashTable* inner_hash = g_hash_table_new_full(
            fq_nmod_hash_func, 
            fq_nmod_equal_func, 
            g_free, // Chaves da tabela interna também serão cópias
            g_array_destroy_wrapper // Valores são GArrays que precisam ser liberados
        );

        for (long j = 0; j < num_polys; j++) {
            // 1. Calcula P_j(x) = y
            fq_nmod_poly_evaluate_fq_nmod(y_eval, polys[j], *x, ctx);

            // 2. Procura pela lista de índices para este 'y'
            GArray* indices = g_hash_table_lookup(inner_hash, y_eval);

            if (indices == NULL) {
                // Se não existe, cria uma nova lista
                indices = g_array_new(FALSE, FALSE, sizeof(long));
                
                // Cria uma cópia de 'y' para usar como chave, pois o 'y_eval' será sobrescrito
                fq_nmod_t* y_key = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
                fq_nmod_init(*y_key, ctx);
                fq_nmod_set(*y_key, y_eval, ctx);
                
                g_hash_table_insert(inner_hash, y_key, indices);
            }
            
            // 3. Adiciona o índice 'j' atual à lista
            g_array_append_val(indices, j);
        }
        
        // Armazena a tabela interna na tabela externa
        // Cria uma cópia de 'x' para a chave
        fq_nmod_t* x_key = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
        fq_nmod_init(*x_key, ctx);
        fq_nmod_set(*x_key, *x, ctx);
        
        g_hash_table_insert(inverted_index, x_key, inner_hash);
    }
    
    fq_nmod_clear(y_eval, ctx);
    return inverted_index;
}

guint fq_nmod_hash_func(gconstpointer key) {
    const fq_nmod_t* element = (const fq_nmod_t*)key;
    return nmod_poly_get_coeff_ui(*element, 0);
}

gboolean fq_nmod_equal_func(gconstpointer a, gconstpointer b) {
    const fq_nmod_t* elem1 = (const fq_nmod_t*)a;
    const fq_nmod_t* elem2 = (const fq_nmod_t*)b;
    return fq_nmod_equal(*elem1, *elem2, global_ctx);
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

void free_matrix(int** matrix, long rows) {
    if (!matrix) return;
    for (long i = 0; i < rows; i++) free(matrix[i]);
    free(matrix);
}

// Função para destruir a GArray dentro da tabela hash interna
void g_array_destroy_wrapper(gpointer data) {
    g_array_free((GArray*)data, TRUE);
}

// Função para destruir a tabela hash interna dentro da tabela externa
void g_hash_table_destroy_wrapper(gpointer data) {
    g_hash_table_destroy((GHashTable*)data);
}