/**
 * @file cff_builder_benchmark.c
 * @brief Implementation of functions for Cover-Free Families (CFFs) construction.
 * 
 * This file contains the functions responsible for generating CFFs using
 * polynomial and monotone constructions over finite fields.
 * 
 * BENCHMARK VERSION: Includes timing measurements for performance analysis.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <time.h>
#include <glib.h>
#include <omp.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_poly.h"
#include "flint/nmod_poly.h"
#include "cff_builder_benchmark.h"
#include "cff_file_generator.h"

/*
 *  GLOBAL VARIABLES
 */

static fq_nmod_ctx_t global_ctx;
static int global_ctx_initialized = 0;

/*
 *  BENCHMARK GLOBAL VARIABLES
 */

/** @brief Flag to enable/disable benchmark mode */
int benchmark_mode = 0;

/** @brief Flag to indicate if this is the last iteration (should save file) */
int benchmark_last_iteration = 0;

/** @brief Accumulated time for Time 1 (inverted index + CFF generation + concatenation) */
double benchmark_time1_accumulated = 0.0;

/** @brief Accumulated time for Time 2 (only CFF matrix generation + concatenation) */
double benchmark_time2_accumulated = 0.0;

/** @brief Accumulated time for generate_single_cff calls */
double benchmark_generate_single_cff_time = 0.0;

/** @brief Accumulated time for concatenation loops */
double benchmark_concatenation_time = 0.0;

/*
 *  FUNCTION PROTOTYPES
 */

/* Main Functions */
void generate_cff(char construction, long fq, long k);
void embed_cff(char construction, int d, long* Fq_steps, long* k_steps);
generated_cffs generate_new_cff_blocks(char construction, int d, long* Fq_steps, long* k_steps, int num_steps);

/* CFF Matrix Generation Functions */
uint64_t** generate_single_cff(long* num_rows, const element_pair* combos, long num_combos, GHashTable* inverted_evals, long num_polys);

/* Hash Table Functions */
guint fq_nmod_hash_func(gconstpointer key);
gboolean fq_nmod_equal_func(gconstpointer a, gconstpointer b);
GHashTable* create_inverted_evaluation_index(long num_points, long num_polys, const fq_nmod_t* points, const fq_nmod_poly_t* polys, const fq_nmod_ctx_t ctx);
void g_array_destroy_wrapper(gpointer data);
void g_hash_table_destroy_wrapper(gpointer data);

/* Element Combination Functions */
combination_partitions generate_combinations(char construction, int dk_size, const subfield_partition* partitions, int num_partitions, const fq_nmod_ctx_t ctx);
void add_pair_to_list(element_pair** list, long* count, long* capacity, const fq_nmod_t x, const fq_nmod_t y, const fq_nmod_ctx_t ctx);

/* Polynomial Functions */
polynomial_partition partition_polynomials(const subfield_partition* partitions, const long* k_steps, int num_steps, const fq_nmod_ctx_t ctx);
fq_nmod_poly_t* generate_polynomials_from_coeffs(long* poly_count, long max_degree, const fq_nmod_t* coeffs, long num_coeffs, const fq_nmod_ctx_t ctx);
void generate_recursive_sorted(fq_nmod_poly_t* poly_list, long* current_index, fq_nmod_poly_t current_poly, long degree, const fq_nmod_t* elements, long num_elements, const fq_nmod_ctx_t ctx);
int fq_nmod_poly_is_in_list(const fq_nmod_poly_t poly, const fq_nmod_poly_t* list, long list_count, const fq_nmod_ctx_t ctx);
void add_poly_to_list(fq_nmod_poly_t** list, long* count, long* capacity, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx);

/* Finite Field Element Functions */
void get_element_by_arithmetic(fq_nmod_t result, ulong i, const fq_nmod_ctx_t ctx);
subfield_partition* partition_by_subfields(const long* Fq_steps, int num_steps, const fq_nmod_ctx_t ctx);
long find_element_index(const fq_nmod_t element, const fq_nmod_t* list, long list_count, const fq_nmod_ctx_t ctx);
void add_element_to_list(fq_nmod_t** list, long* count, long* capacity, const fq_nmod_t element, const fq_nmod_ctx_t ctx);

/* Mathematical Utility Functions */
static int is_prime(long n);
static int decompose_prime_power(long q, long* p_out, long* n_out);

/* Memory Deallocation Functions */
void free_matrix(uint64_t** matrix, long rows);
void free_subfield_partitions(subfield_partition* partitions, int num_steps, const fq_nmod_ctx_t ctx);
void free_combination_partitions(combination_partitions* combos, const fq_nmod_ctx_t ctx);
static void free_generated_cffs(generated_cffs* cffs);
static void free_polynomial_partition(polynomial_partition* poly_part, const fq_nmod_ctx_t ctx);

/* 
 *  BENCHMARK HELPER FUNCTIONS
 */

/**
 * @brief Gets current time in seconds with high precision.
 * @return Current time in seconds.
 */
static inline double get_time_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

/**
 * @brief Resets benchmark timing accumulators for a new iteration.
 */
void benchmark_reset_iteration(void) {
    benchmark_generate_single_cff_time = 0.0;
    benchmark_concatenation_time = 0.0;
}

/* 
 *  MAIN FUNCTIONS
 */

/**
 * @brief Generates an initial CFF from basic parameters.
 * 
 * Creates a CFF from scratch using the provided parameters and saves it to file.
 * 
 * @param construction Construction type ('p' for polynomial, 'm' for monotone).
 * @param fq Finite field size.
 * @param k Maximum polynomial degree.
 */
void generate_cff(char construction, long fq, long k) {
    double time1_start = 0.0, time1_end = 0.0;
    
    /* Reset iteration timers */
    if (benchmark_mode) {
        benchmark_reset_iteration();
    }
    
    int d;
    long t, n;
    
    d = (fq-1)/k;
    t = fq*fq;
    n = (long)pow(fq, k + 1);
    
    char filename0[100];
    snprintf(filename0, sizeof(filename0), "CFFs/%d-CFF(%ld,%ld).txt", d, t, n);

    long fq_array[1] = { fq };
    long k_array[1] = { k };
    int num_steps = 1;

    /* TIME 1 START: Inverted index + CFF generation */
    if (benchmark_mode) {
        time1_start = get_time_seconds();
    }
    
    generated_cffs new_blocks = generate_new_cff_blocks(construction, d, fq_array, k_array, num_steps);
    
    /* TIME 1 END */
    if (benchmark_mode) {
        time1_end = get_time_seconds();
        benchmark_time1_accumulated += (time1_end - time1_start);
        /* For generate_cff, Time 2 = time spent in generate_single_cff */
        benchmark_time2_accumulated += benchmark_generate_single_cff_time;
    }

    uint64_t** final_cff = new_blocks.cff_new;
    long final_rows = new_blocks.rows_new;
    long final_cols = new_blocks.cols_new;

    if (final_cff == NULL) {
        printf("Error: Failed to generate initial CFF matrix.\n");
        return;
    }

    /* Only write to file on last iteration or if not in benchmark mode */
    if (!benchmark_mode || benchmark_last_iteration) {
        write_cff_to_file(filename0, construction, d, fq_array, 1, k_array, 1, final_cff, final_rows, final_cols);
    }

    new_blocks.cff_old_new = NULL;
    new_blocks.cff_new_old = NULL;

    free_generated_cffs(&new_blocks); 
}

/**
 * @brief Performs embedding of an existing CFF to a larger field.
 * 
 * Reads an existing CFF from file and expands it to a larger finite field,
 * generating the necessary new blocks and saving the expanded CFF.
 * 
 * @param construction Construction type ('p' for polynomial, 'm' for monotone).
 * @param d CFF parameter d.
 * @param Fq_steps Array with finite field sizes.
 * @param k_steps Array with maximum polynomial degrees.
 */
void embed_cff(char construction, int d, long* Fq_steps, long* k_steps){
    double time1_start = 0.0, time1_end = 0.0;
    double concat_start = 0.0, concat_end = 0.0;
    
    /* Reset iteration timers */
    if (benchmark_mode) {
        benchmark_reset_iteration();
    }
    
    int d0;
    long t0, n0;
    
    if (construction == 'm') {
        d0 = d;
        t0 = (d * k_steps[0] + 1) * Fq_steps[0];
        n0 = (long)pow(Fq_steps[0], k_steps[0] + 1);
    } else {
        d0 = (Fq_steps[0]-1)/k_steps[0];
        t0 = Fq_steps[0] * Fq_steps[0];
        n0 = (long)pow(Fq_steps[0], k_steps[0] + 1);
    }
    
    char filename0[100]; 
    snprintf(filename0, sizeof(filename0), "CFFs/%d-CFF(%ld,%ld).txt", d0, t0, n0);

    long old_rows = 0, old_cols = 0;
    struct cff_parameters* params = read_parameters(filename0);
    if (params == NULL) {
        printf("Error reading parameters from file %s\n", filename0);
        return;
    }

    uint64_t** cff_old_old = read_cff_from_file(filename0, &old_rows, &old_cols);

    int new_fqs_count = params->fqs_count + 1;
    int new_ks_count = params->ks_count + 1;
    long* new_Fq_steps = (long*) malloc(new_fqs_count * sizeof(long));
    long* new_k_steps  = (long*) malloc(new_ks_count * sizeof(long));

    if (new_Fq_steps == NULL || new_k_steps == NULL) {
        printf("Error: Failed to allocate memory for step conversion.\n");
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

    /* TIME 1 START: Inverted index + CFF generation + concatenation */
    if (benchmark_mode) {
        time1_start = get_time_seconds();
    }

    generated_cffs new_blocks = generate_new_cff_blocks(construction, d0, new_Fq_steps, new_k_steps, new_fqs_count);

    long new_total_rows = old_rows + new_blocks.rows_new_old;

    long width_top = old_cols + new_blocks.cols_old_new;
    long width_bottom = new_blocks.cols_new_old + new_blocks.cols_new;
    long new_total_cols = (width_top > width_bottom) ? width_top : width_bottom;
    
    long words_per_row = WORDS_FOR_BITS(new_total_cols);
    uint64_t** final_cff = (uint64_t**) malloc(new_total_rows * sizeof(uint64_t*));
    for (long i = 0; i < new_total_rows; i++) {
        final_cff[i] = (uint64_t*) calloc(words_per_row, sizeof(uint64_t));
    }

    /* CONCATENATION START - also measured for Time 2 */
    if (benchmark_mode) {
        concat_start = get_time_seconds();
    }

    for(long i=0; i < old_rows; i++) {
        for (long j = 0; j < old_cols; j++) {
            if (GET_BIT(cff_old_old[i], j)) {
                SET_BIT(final_cff[i], j);
            }
        }
    }

    for(long i=0; i < new_blocks.rows_old_new; i++) {
        for (long j = 0; j < new_blocks.cols_old_new; j++) {
            if (GET_BIT(new_blocks.cff_old_new[i], j)) {
                SET_BIT(final_cff[i], old_cols + j);
            }
        }
    }
    
    for(long i=0; i < new_blocks.rows_new_old; i++) {
        for (long j = 0; j < new_blocks.cols_new_old; j++) {
            if (GET_BIT(new_blocks.cff_new_old[i], j)) {
                SET_BIT(final_cff[old_rows + i], j);
            }
        }
    }
    
    for(long i=0; i < new_blocks.rows_new; i++) {
        for (long j = 0; j < new_blocks.cols_new; j++) {
            if (GET_BIT(new_blocks.cff_new[i], j)) {
                SET_BIT(final_cff[old_rows + i], new_blocks.cols_new_old + j);
            }
        }
    }

    /* CONCATENATION END */
    if (benchmark_mode) {
        concat_end = get_time_seconds();
        benchmark_concatenation_time += (concat_end - concat_start);
    }

    /* TIME 1 END */
    if (benchmark_mode) {
        time1_end = get_time_seconds();
        benchmark_time1_accumulated += (time1_end - time1_start);
        /* Time 2 = generate_single_cff time + concatenation time */
        benchmark_time2_accumulated += benchmark_generate_single_cff_time + benchmark_concatenation_time;
    }

    int d1;
    long t1, n1;
    
    if (construction == 'm') {
        d1 = d;
        t1 = (d * k_steps[0] + 1) * Fq_steps[1];
        n1 = (long)pow(Fq_steps[1], k_steps[0] + 1);
    } else {
        d1 = (Fq_steps[1]-1)/k_steps[1];
        t1 = Fq_steps[1] * Fq_steps[1];
        n1 = (long)pow(Fq_steps[1], k_steps[1] + 1);
    }
    
    char filename1[100]; 
    snprintf(filename1, sizeof(filename1), "CFFs/%d-CFF(%ld,%ld).txt", d1, t1, n1);

    /* Only write to file on last iteration or if not in benchmark mode */
    if (!benchmark_mode || benchmark_last_iteration) {
        write_cff_to_file(filename1, construction, d1, new_Fq_steps, new_fqs_count, new_k_steps, new_ks_count, final_cff, new_total_rows, new_total_cols);
    }

    free_matrix(cff_old_old, old_rows);
    free_generated_cffs(&new_blocks);
    free_matrix(final_cff, new_total_rows);
    free(new_Fq_steps);
    free(new_k_steps);
    free(params->Fqs);
    free(params->ks);
    free(params);
}

/**
 * @brief Generates new CFF blocks for an embedding step.
 * 
 * Main function that orchestrates the generation of the three blocks needed
 * to expand a CFF: old_new, new_old, and new_new.
 * 
 * @param construction Construction type ('p' or 'm').
 * @param d CFF parameter d (used for monotone construction).
 * @param Fq_steps Array with finite field sizes.
 * @param k_steps Array with maximum polynomial degrees.
 * @param num_steps Number of steps.
 * @return Structure containing the three generated blocks.
 */
generated_cffs generate_new_cff_blocks(char construction, int d, long* Fq_steps, long* k_steps, int num_steps) {
    generated_cffs result = {0};
    
    long q_final = Fq_steps[num_steps - 1];  
    long p, n;
    
    if (!decompose_prime_power(q_final, &p, &n)) {
        fprintf(stderr, "Error: %ld is not a prime power!\n", q_final);
        return result;
    }

    fmpz_t pz;
    fmpz_init(pz);
    fmpz_set_ui(pz, (ulong)p);

    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx, (ulong)p, (slong)n, "a");
    
    subfield_partition* partitions = partition_by_subfields(Fq_steps, num_steps, ctx);

    polynomial_partition poly_part = partition_polynomials(partitions, k_steps, num_steps, ctx);

    long num_new_rows = 0;
    int dk_size = 0;
    if(construction == 'p'){
        if (num_steps == 1) {
            num_new_rows = (d*k_steps[0]+1)*Fq_steps[0];
        } else {
            num_new_rows = (Fq_steps[num_steps-1] * Fq_steps[num_steps-1]) - (Fq_steps[num_steps-2] * Fq_steps[num_steps-2]);
        }
        dk_size = Fq_steps[num_steps-1];
    } else if (construction == 'm') {
        num_new_rows = (d*k_steps[0]+1)*Fq_steps[num_steps-1] - (d*k_steps[num_steps-2]+1)*Fq_steps[num_steps-2];
        dk_size = d*k_steps[0]+1;
    };  

    combination_partitions combos = generate_combinations(construction, dk_size, partitions, num_steps, ctx);

    subfield_partition all_partition = partitions[num_steps - 1];
    fq_nmod_t* points_for_eval = all_partition.all_elements;
    long num_points = all_partition.count_all;
    GHashTable* inverted_index_old = create_inverted_evaluation_index(num_points, poly_part.num_old_polys, points_for_eval, poly_part.old_polys, ctx);
    GHashTable* inverted_index_new = create_inverted_evaluation_index(num_points, poly_part.num_new_polys, points_for_eval, poly_part.new_polys, ctx);

    if (!global_ctx_initialized) {
        fq_nmod_ctx_init_modulus(global_ctx, fq_nmod_ctx_modulus(ctx), "a");
        global_ctx_initialized = 1;
    }
    
    /* Timing for generate_single_cff calls */
    double gen_start, gen_end;
    
    if (benchmark_mode) {
        gen_start = get_time_seconds();
    }
    result.cff_old_new = generate_single_cff(&result.rows_old_new, combos.combos_old, combos.count_old, inverted_index_new, poly_part.num_new_polys);
    if (benchmark_mode) {
        gen_end = get_time_seconds();
        benchmark_generate_single_cff_time += (gen_end - gen_start);
    }
    result.cols_old_new = poly_part.num_new_polys;
    
    if (benchmark_mode) {
        gen_start = get_time_seconds();
    }
    result.cff_new_old = generate_single_cff(&result.rows_new_old, combos.combos_new, num_new_rows, inverted_index_old, poly_part.num_old_polys);
    if (benchmark_mode) {
        gen_end = get_time_seconds();
        benchmark_generate_single_cff_time += (gen_end - gen_start);
    }
    result.cols_new_old = poly_part.num_old_polys;

    if (benchmark_mode) {
        gen_start = get_time_seconds();
    }
    result.cff_new = generate_single_cff(&result.rows_new, combos.combos_new, num_new_rows, inverted_index_new, poly_part.num_new_polys);
    if (benchmark_mode) {
        gen_end = get_time_seconds();
        benchmark_generate_single_cff_time += (gen_end - gen_start);
    }
    result.cols_new = poly_part.num_new_polys;
    
    g_hash_table_destroy(inverted_index_old);
    g_hash_table_destroy(inverted_index_new);
    free_combination_partitions(&combos, ctx);
    free_polynomial_partition(&poly_part, ctx);
    free_subfield_partitions(partitions, num_steps, ctx);
    fq_nmod_ctx_clear(ctx);
    fq_nmod_ctx_clear(global_ctx);
    global_ctx_initialized = 0;

    return result;
}

/*
 *  CFF MATRIX GENERATION FUNCTION
 */

/**
 * @brief Generates a CFF matrix (bitmap) from combinations and evaluations.
 * 
 * Creates a binary matrix where each row corresponds to a pair (x, y) and
 * each column corresponds to a polynomial. Bit (i, j) is 1 if polynomial j
 * evaluates to y at point x of pair i.
 * 
 * @param num_rows Pointer to store the number of rows.
 * @param combos Array of element pairs.
 * @param num_combos Number of pairs.
 * @param inverted_evals Inverted evaluation index.
 * @param num_polys Number of polynomials (columns).
 * @return CFF matrix in bitmap format.
 */
uint64_t** generate_single_cff(long* num_rows, const element_pair* combos, long num_combos, GHashTable* inverted_evals, long num_polys) {
    *num_rows = num_combos;
    if (num_combos == 0) return NULL;

    long words_per_row = WORDS_FOR_BITS(num_polys);

    uint64_t** cff_matrix = (uint64_t**) malloc(num_combos * sizeof(uint64_t*));
    if (cff_matrix == NULL) exit(EXIT_FAILURE);

    #pragma omp parallel for schedule(dynamic)
    for (long i = 0; i < num_combos; i++) {

        const fq_nmod_t* x = &combos[i].x;
        const fq_nmod_t* y = &combos[i].y;

        cff_matrix[i] = (uint64_t*) calloc(words_per_row, sizeof(uint64_t));

        GHashTable* inner_hash = g_hash_table_lookup(inverted_evals, x);

        if (inner_hash != NULL) {
            GArray* indices = g_hash_table_lookup(inner_hash, y);

            if (indices != NULL) {
                for (guint k = 0; k < indices->len; k++) {
                    long j = g_array_index(indices, long, k);
                    SET_BIT(cff_matrix[i], j);
                }
            }
        }
    }

    return cff_matrix;
}

/* 
 *  HASH TABLE FUNCTIONS
 */

/**
 * @brief Creates an inverted index of polynomial evaluations.
 * 
 * For each point x, creates a mapping y -> list of polynomial indices
 * where p(x) = y. This enables fast queries of which polynomials evaluate
 * to a given pair (x, y).
 * 
 * @param num_points Number of evaluation points.
 * @param num_polys Number of polynomials.
 * @param points Array of evaluation points.
 * @param polys Array of polynomials.
 * @param ctx Finite field context.
 * @return Hash table with the inverted index.
 */
GHashTable* create_inverted_evaluation_index(long num_points, long num_polys, const fq_nmod_t* points, const fq_nmod_poly_t* polys, const fq_nmod_ctx_t ctx) {
    GHashTable* inverted_index = g_hash_table_new_full(fq_nmod_hash_func, fq_nmod_equal_func, g_free, g_hash_table_destroy_wrapper);

    #pragma omp parallel
    {
        fq_nmod_t y_eval_local;
        fq_nmod_init(y_eval_local, ctx);

        fq_nmod_t* x_key_copy; 
        fq_nmod_t* y_key_copy;

        #pragma omp for schedule(dynamic)
        for (long i = 0; i < num_points; i++) {
            const fq_nmod_t* x = &points[i];

            GHashTable* inner_hash_local = g_hash_table_new_full(fq_nmod_hash_func, fq_nmod_equal_func, g_free, g_array_destroy_wrapper);

            for (long j = 0; j < num_polys; j++) {
                fq_nmod_poly_evaluate_fq_nmod(y_eval_local, polys[j], *x, ctx);

                GArray* indices = g_hash_table_lookup(inner_hash_local, y_eval_local);
                if (indices == NULL) {
                    indices = g_array_new(FALSE, FALSE, sizeof(long));
                    y_key_copy = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
                    fq_nmod_init(*y_key_copy, ctx);
                    fq_nmod_set(*y_key_copy, y_eval_local, ctx);
                    g_hash_table_insert(inner_hash_local, y_key_copy, indices);
                }
                g_array_append_val(indices, j);
            }

            x_key_copy = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
            fq_nmod_init(*x_key_copy, ctx);
            fq_nmod_set(*x_key_copy, *x, ctx);

            #pragma omp critical
            g_hash_table_insert(inverted_index, x_key_copy, inner_hash_local);
        }
        fq_nmod_clear(y_eval_local, ctx);
    }

    return inverted_index;
}

/**
 * @brief Hash function for fq_nmod_t elements.
 * 
 * Computes a hash value based on the degree-zero coefficient of the element.
 * 
 * @param key Pointer to the fq_nmod_t element.
 * @return Computed hash value.
 */
guint fq_nmod_hash_func(gconstpointer key) {
    const fq_nmod_t* element = (const fq_nmod_t*)key;
    return nmod_poly_get_coeff_ui(*element, 0);
}

/**
 * @brief Equality comparison function for fq_nmod_t elements.
 * 
 * Compares two finite field elements using the global context.
 * 
 * @param a Pointer to the first element.
 * @param b Pointer to the second element.
 * @return TRUE if elements are equal, FALSE otherwise.
 */
gboolean fq_nmod_equal_func(gconstpointer a, gconstpointer b) {
    const fq_nmod_t* elem1 = (const fq_nmod_t*)a;
    const fq_nmod_t* elem2 = (const fq_nmod_t*)b;
    return fq_nmod_equal(*elem1, *elem2, global_ctx);
}

/**
 * @brief Wrapper to destroy a GArray inside a hash table.
 * 
 * @param data Pointer to the GArray to be destroyed.
 */
void g_array_destroy_wrapper(gpointer data) {
    g_array_free((GArray*)data, TRUE);
}

/**
 * @brief Wrapper to destroy an inner hash table.
 * 
 * @param data Pointer to the GHashTable to be destroyed.
 */
void g_hash_table_destroy_wrapper(gpointer data) {
    g_hash_table_destroy((GHashTable*)data);
}

/* 
 *  COMBINATION HELPER FUNCTIONS
 */

/**
 * @brief Generates lists of "old" and "new" element pairs.
 * 
 * Creates combinations of pairs (x, y) based on subfield partitions,
 * classifying them as old (from previous steps) or new.
 * 
 * @param construction Construction type ('p' for polynomial, 'm' for monotone).
 * @param dk_size Size of dk block (used only for monotone construction).
 * @param partitions Array of subfield partitions.
 * @param num_partitions Number of partitions.
 * @param ctx Finite field context.
 * @return Structure containing the partitioned combinations.
 */
combination_partitions generate_combinations(char construction, int dk_size, const subfield_partition* partitions, int num_partitions, const fq_nmod_ctx_t ctx) {
    combination_partitions result = {0};

    fq_nmod_t* all_accumulated_elements = NULL;
    long accumulated_count = 0;
    long accumulated_capacity = 0;
    fq_nmod_t* dk_block_elements = NULL; 

    if(construction == 'm'){
        dk_block_elements = malloc(dk_size * sizeof(fq_nmod_t));
        int start_index = 0; 
        
        fq_nmod_t* current_only_elements = partitions[0].only_elements;
        for(int i=0; i<dk_size; i++){
            fq_nmod_init(dk_block_elements[i], ctx);
            fq_nmod_set(dk_block_elements[i], current_only_elements[start_index + i], ctx);
        }
    }

    for (int i = 0; i < num_partitions; i++) {
        fq_nmod_t* current_only_elements = partitions[i].only_elements;
        long current_only_count = partitions[i].count_only;

        element_pair** target_list = (i == num_partitions - 1) ? &result.combos_new : &result.combos_old;
        long* target_count = (i == num_partitions - 1) ? &result.count_new : &result.count_old;
        long* target_capacity = (i == num_partitions - 1) ? &result.capacity_new : &result.capacity_old;

        if (i == 0) {
            for (long ix = 0; ix < current_only_count; ix++) {
                for (long iy = 0; iy < current_only_count; iy++) {
                    add_pair_to_list(target_list, target_count, target_capacity, current_only_elements[ix], current_only_elements[iy], ctx);
                }
            }
        } else {
            if(construction == 'p'){
                for (long ix = 0; ix < accumulated_count; ix++) {
                    for (long iy = 0; iy < current_only_count; iy++) {
                        add_pair_to_list(target_list, target_count, target_capacity, all_accumulated_elements[ix], current_only_elements[iy], ctx);
                    }
                }
                for (long j = 0; j < current_only_count; j++) {
                    add_element_to_list(&all_accumulated_elements, &accumulated_count, &accumulated_capacity, current_only_elements[j], ctx);
                }
                for (long ix = 0; ix < current_only_count; ix++) {
                    for (long iy = 0; iy < accumulated_count; iy++) {
                        add_pair_to_list(target_list, target_count, target_capacity, current_only_elements[ix], all_accumulated_elements[iy], ctx);
                    }
                } 
            } else if (construction == 'm'){
                for (long ix = 0; ix < dk_size; ix++) {
                    for (long iy = 0; iy < current_only_count; iy++) {
                        add_pair_to_list(target_list, target_count, target_capacity, dk_block_elements[ix], current_only_elements[iy], ctx);
                    }
                }
            }
        }
        
        if (i == 0) {
             for (long j = 0; j < current_only_count; j++) {
                add_element_to_list(&all_accumulated_elements, &accumulated_count, &accumulated_capacity, current_only_elements[j], ctx);
            }
        }
    }
    
    for (long i = 0; i < accumulated_count; i++) {
        fq_nmod_clear(all_accumulated_elements[i], ctx);
    }
    free(all_accumulated_elements);

    return result;
}

/**
 * @brief Adds an element pair (x, y) to a dynamic array of pairs.
 * 
 * Automatically manages memory reallocation when needed,
 * doubling the array capacity.
 * 
 * @param list Pointer to the pair array.
 * @param count Pointer to the pair counter.
 * @param capacity Pointer to the current array capacity.
 * @param x First element of the pair.
 * @param y Second element of the pair.
 * @param ctx Finite field context.
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

/* =============================================================================
 * POLYNOMIAL HELPER FUNCTIONS
 * ========================================================================== */

/**
 * @brief Partitions polynomials into "old" and "new" based on subfields.
 * 
 * Generates polynomials for each subfield step and classifies them as old
 * (belonging to previous steps) or new (exclusive to the last step).
 * 
 * @param partitions Array of subfield partitions.
 * @param k_steps Array with maximum degrees for each step.
 * @param num_steps Number of steps.
 * @param ctx Finite field context.
 * @return Structure containing the partitioned polynomials.
 */
polynomial_partition partition_polynomials(const subfield_partition* partitions, const long* k_steps, int num_steps, const fq_nmod_ctx_t ctx) {
    polynomial_partition result = {0};
    
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
                    fprintf(stderr, "Error reallocating memory!\n");
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

    result.old_polys = old_polys;
    result.num_old_polys = num_old_polys_total;
    result.new_polys = new_polys;
    result.num_new_polys = num_new_polys;
    result.all_polys = all_polys;
    result.num_all_polys = num_all_polys;

    return result;
}

/**
 * @brief Generates all polynomials up to a maximum degree with given coefficients.
 * 
 * @param poly_count Pointer to store the number of generated polynomials.
 * @param max_degree Maximum polynomial degree.
 * @param coeffs Array of possible coefficients.
 * @param num_coeffs Number of possible coefficients.
 * @param ctx Finite field context.
 * @return Array of generated polynomials.
 */
fq_nmod_poly_t* generate_polynomials_from_coeffs(long* poly_count, long max_degree, const fq_nmod_t* coeffs, long num_coeffs, const fq_nmod_ctx_t ctx) {
    fmpz_t num_coeffs_z, total_polys_z;
    fmpz_init(num_coeffs_z);
    fmpz_init(total_polys_z);
    fmpz_set_si(num_coeffs_z, num_coeffs);

    fmpz_pow_ui(total_polys_z, num_coeffs_z, max_degree + 1);
    *poly_count = fmpz_get_si(total_polys_z);

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
 * @brief Recursive function that builds polynomials prioritizing higher degree terms.
 * 
 * Generates polynomials in an ordered manner, starting from the highest degree
 * coefficient and decreasing to the constant term.
 * 
 * @param poly_list Array to store the generated polynomials.
 * @param current_index Pointer to the current index in the list.
 * @param current_poly Polynomial being built.
 * @param degree Current degree being processed.
 * @param elements Array of possible elements as coefficients.
 * @param num_elements Number of possible elements.
 * @param ctx Finite field context.
 */
void generate_recursive_sorted(
    fq_nmod_poly_t* poly_list, 
    long* current_index, 
    fq_nmod_poly_t current_poly, 
    long degree,
    const fq_nmod_t* elements, 
    long num_elements,
    const fq_nmod_ctx_t ctx) 
{
    if (degree < 0) {
        fq_nmod_poly_set(poly_list[*current_index], current_poly, ctx);
        (*current_index)++;
        return;
    }

    for (long i = 0; i < num_elements; i++) {
        fq_nmod_poly_set_coeff(current_poly, degree, elements[i], ctx);
        generate_recursive_sorted(poly_list, current_index, current_poly, degree - 1, elements, num_elements, ctx);
    }
}

/**
 * @brief Checks if a polynomial exists in an array of polynomials.
 * 
 * @param poly Polynomial to search for.
 * @param list Array of polynomials.
 * @param list_count Number of polynomials in the array.
 * @param ctx Finite field context.
 * @return 1 if found, 0 otherwise.
 */
int fq_nmod_poly_is_in_list(const fq_nmod_poly_t poly, const fq_nmod_poly_t* list, long list_count, const fq_nmod_ctx_t ctx) {
    for (long i = 0; i < list_count; i++) {
        if (fq_nmod_poly_equal(poly, list[i], ctx)) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Adds a fq_nmod_poly_t polynomial to a dynamic array.
 * 
 * Automatically manages memory reallocation when needed,
 * doubling the array capacity.
 * 
 * @param list Pointer to the polynomial array.
 * @param count Pointer to the polynomial counter.
 * @param capacity Pointer to the current array capacity.
 * @param poly Polynomial to be added.
 * @param ctx Finite field context.
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

/* 
 *  FINITE FIELD ELEMENT HELPER FUNCTIONS
 */

/**
 * @brief Constructs a finite field element from its arithmetic index.
 * 
 * Converts an integer index to its representation as a finite field element
 * using the characteristic prime p as base.
 * 
 * @param result Resulting element.
 * @param i Element index (0 to q-1).
 * @param ctx Finite field context.
 */
void get_element_by_arithmetic(fq_nmod_t result, ulong i, const fq_nmod_ctx_t ctx) {
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
 * @brief Partitions finite field elements by subfield membership.
 * 
 * Iterates through all elements of the main field and classifies them
 * according to their membership in the subfields specified in Fq_steps.
 * 
 * @param Fq_steps Array with subfield sizes.
 * @param num_steps Number of subfields.
 * @param ctx Finite field context.
 * @return Array of subfield partitions.
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
 * @brief Finds the index of an element in an element array.
 * 
 * @param element Element to search for.
 * @param list Array of elements.
 * @param list_count Number of elements in the array.
 * @param ctx Finite field context.
 * @return Element index if found, -1 otherwise.
 */
long find_element_index(const fq_nmod_t element, const fq_nmod_t* list, long list_count, const fq_nmod_ctx_t ctx) {
    for (long i = 0; i < list_count; i++) {
        if (fq_nmod_equal(element, list[i], ctx)) {
            return i;
        }
    }
    return -1;
}

/**
 * @brief Adds a fq_nmod_t element to a dynamic array.
 * 
 * Automatically manages memory reallocation when needed,
 * doubling the array capacity.
 * 
 * @param list Pointer to the element array.
 * @param count Pointer to the element counter.
 * @param capacity Pointer to the current array capacity.
 * @param element Element to be added.
 * @param ctx Finite field context.
 */
void add_element_to_list(fq_nmod_t** list, long* count, long* capacity, const fq_nmod_t element, const fq_nmod_ctx_t ctx) {
    if (*count >= *capacity) {
        *capacity = (*capacity == 0) ? 8 : (*capacity) * 2;
        *list = realloc(*list, (*capacity) * sizeof(fq_nmod_t));
        for (long i = *count; i < *capacity; i++) {
            fq_nmod_init((*list)[i], ctx);
        }
    }
    fq_nmod_set((*list)[*count], element, ctx);
    (*count)++;
}

/*
 * MATHEMATICAL UTILITY FUNCTIONS
 */

/**
 * @brief Checks if a number is prime.
 * 
 * @param n Number to check.
 * @return 1 if n is prime, 0 otherwise.
 */
static int is_prime(long n) {
    if (n <= 1) return 0;
    if (n <= 3) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (long i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return 0;
    }
    return 1;
}

/**
 * @brief Decomposes q into p^n where p is prime.
 * 
 * @param q Finite field size (must be a prime power).
 * @param p_out Pointer to store the prime characteristic.
 * @param n_out Pointer to store the exponent.
 * @return 1 on success, 0 if q is not a prime power.
 */
static int decompose_prime_power(long q, long* p_out, long* n_out) {
    if (q <= 1) return 0;
    
    if (is_prime(q)) {
        *p_out = q;
        *n_out = 1;
        return 1;
    }
    
    for (long p = 2; p * p <= q; p++) {
        if (!is_prime(p)) continue;
        
        long temp = q;
        long n = 0;
        
        while (temp % p == 0) {
            temp /= p;
            n++;
        }
        
        if (temp == 1 && n > 0) {
            *p_out = p;
            *n_out = n;
            return 1;
        }
    }
    
    return 0;
}

/* 
 *  MEMORY DEALLOCATION FUNCTIONS
 */

/**
 * @brief Frees the memory of a uint64_t matrix.
 * 
 * @param matrix Matrix to be freed.
 * @param rows Number of rows in the matrix.
 */
void free_matrix(uint64_t** matrix, long rows) {
    if (!matrix) return;
    for (long i = 0; i < rows; i++) free(matrix[i]);
    free(matrix);
}

/**
 * @brief Frees the memory of a subfield partition array.
 * 
 * @param partitions Array of partitions to be freed.
 * @param num_steps Number of partitions in the array.
 * @param ctx Finite field context.
 */
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

/**
 * @brief Frees the memory of a combination partition structure.
 * 
 * @param combos Pointer to the structure to be freed.
 * @param ctx Finite field context.
 */
void free_combination_partitions(combination_partitions* combos, const fq_nmod_ctx_t ctx) {
    if (!combos) return;
    for (long i = 0; i < combos->capacity_old; i++) { fq_nmod_clear(combos->combos_old[i].x, ctx); fq_nmod_clear(combos->combos_old[i].y, ctx); }
    free(combos->combos_old);
    for (long i = 0; i < combos->capacity_new; i++) { fq_nmod_clear(combos->combos_new[i].x, ctx); fq_nmod_clear(combos->combos_new[i].y, ctx); }
    free(combos->combos_new);
}

/**
 * @brief Frees the memory of a generated CFFs structure.
 * 
 * @param cffs Pointer to the structure to be freed.
 */
static void free_generated_cffs(generated_cffs* cffs) {
    if (!cffs) return;
    if (cffs->cff_old_new) { for (long i = 0; i < cffs->rows_old_new; i++) free(cffs->cff_old_new[i]); free(cffs->cff_old_new); }
    if (cffs->cff_new_old) { for (long i = 0; i < cffs->rows_new_old; i++) free(cffs->cff_new_old[i]); free(cffs->cff_new_old); }
    if (cffs->cff_new) { for (long i = 0; i < cffs->rows_new; i++) free(cffs->cff_new[i]); free(cffs->cff_new); }
}

/**
 * @brief Frees the memory of a polynomial partition structure.
 * 
 * @param poly_part Pointer to the structure to be freed.
 * @param ctx Finite field context.
 */
static void free_polynomial_partition(polynomial_partition* poly_part, const fq_nmod_ctx_t ctx) {
    if (!poly_part) return;
    if (poly_part->old_polys != NULL) { 
        for (long i = 0; i < poly_part->num_old_polys; i++) fq_nmod_poly_clear(poly_part->old_polys[i], ctx); 
        free(poly_part->old_polys); 
    }
    if (poly_part->new_polys != NULL) { 
        for (long i = 0; i < poly_part->num_new_polys; i++) fq_nmod_poly_clear(poly_part->new_polys[i], ctx); 
        free(poly_part->new_polys); 
    }
    if (poly_part->all_polys != NULL) { 
        for (long i = 0; i < poly_part->num_all_polys; i++) fq_nmod_poly_clear(poly_part->all_polys[i], ctx); 
        free(poly_part->all_polys); 
    }
}