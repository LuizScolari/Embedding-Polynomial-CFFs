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

// Struct para um par de elementos (x, y)
typedef struct {
    fq_nmod_t x;
    fq_nmod_t y;
} element_pair;

// Struct para conter as duas listas de combinações resultantes
typedef struct {
    element_pair* combos_old;
    long count_old;
    long capacity_old;
    element_pair* combos_new;
    long count_new;
    long capacity_new;
} combination_partitions;


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
void free_cff_matrix(int** matrix, long num_rows);
void print_cff_matrix(int** matrix, long num_rows, long num_cols);


int main() {
    // --- ETAPA 0: PARÂMETROS E INICIALIZAÇÃO ---
    fmpz_t pz; fmpz_init(pz); fmpz_set_ui(pz, 2);
    long D = 2; // F_16
    long Fq_steps[] = {2, 4}; 
    int num_steps = 2;

    nmod_poly_t mod_poly; nmod_poly_init(mod_poly, 2);
    nmod_poly_set_coeff_ui(mod_poly, 2, 1); 
    nmod_poly_set_coeff_ui(mod_poly, 1, 1); 
    nmod_poly_set_coeff_ui(mod_poly, 0, 1);
    fq_nmod_ctx_t ctx; fq_nmod_ctx_init_modulus(ctx, mod_poly, "a");
    nmod_poly_clear(mod_poly);
    
    // --- ETAPA 1: PARTICIONAR ELEMENTOS ---
    subfield_partition* partitions = partition_by_subfields(Fq_steps, num_steps, ctx);

    // --- ETAPA 2: PARTICIONAR POLINÔMIOS ('old' vs 'new') ---
    long max_poly_degree = 1;
    // 'old_polys' são os polinômios do penúltimo subcorpo (F_4)
    subfield_partition old_partition = partitions[num_steps - 2];
    long num_old_polys = 0;
    fq_nmod_poly_t* old_polys = generate_polynomials_from_coeffs(&num_old_polys, max_poly_degree, old_partition.all_elements, old_partition.count_all, ctx);
    
    // 'all_polys' são os polinômios do corpo inteiro (F_16)
    subfield_partition all_partition = partitions[num_steps - 1];
    long num_all_polys = 0;
    fq_nmod_poly_t* all_polys = generate_polynomials_from_coeffs(&num_all_polys, max_poly_degree, all_partition.all_elements, all_partition.count_all, ctx);
    
    // 'new_polys' são os que estão em 'all_polys' mas não em 'old_polys'
    fq_nmod_poly_t* new_polys = (fq_nmod_poly_t*) malloc(num_all_polys * sizeof(fq_nmod_poly_t));
    long num_new_polys = 0;
    for (long i = 0; i < num_all_polys; i++) {
        if (!fq_nmod_poly_is_in_list(all_polys[i], old_polys, num_old_polys, ctx)) {
            fq_nmod_poly_init(new_polys[num_new_polys], ctx);
            fq_nmod_poly_set(new_polys[num_new_polys], all_polys[i], ctx);
            num_new_polys++;
        }
    }

    // --- ETAPA 3: GERAR COMBINAÇÕES DE ELEMENTOS ---
    combination_partitions combos = generate_combinations(partitions, num_steps, ctx);

    // --- ETAPA 4: CRIAR MATRIZES DE AVALIAÇÃO ---
    fq_nmod_t* points_for_eval = all_partition.all_elements;
    long num_points = all_partition.count_all;

    fq_nmod_t** evals_old = create_evaluation_matrix(num_points, num_old_polys, points_for_eval, old_polys, ctx);
    fq_nmod_t** evals_new = create_evaluation_matrix(num_points, num_new_polys, points_for_eval, new_polys, ctx);

    // --- ETAPA 5: GERAR MATRIZES CFF FINAIS ---
    long rows_cff1 = 0;
    int** cff_old_new = generate_single_cff(&rows_cff1, combos.combos_old, combos.count_old, points_for_eval, num_points, evals_new, num_new_polys, ctx);
    
    long rows_cff2 = 0;
    int** cff_new_old = generate_single_cff(&rows_cff2, combos.combos_new, combos.count_new, points_for_eval, num_points, evals_old, num_old_polys, ctx);

    long rows_cff3 = 0;
    int** cff_new = generate_single_cff(&rows_cff3, combos.combos_new, combos.count_new, points_for_eval, num_points, evals_new, num_new_polys, ctx);

    // --- ETAPA 6: IMPRIMIR RESULTADOS ---
    print_cff_matrix(cff_old_new, rows_cff1, num_new_polys);
    print_cff_matrix(cff_new_old, rows_cff2, num_old_polys);
    print_cff_matrix(cff_new, rows_cff3, num_new_polys);
    
    // --- ETAPA 7: LIMPEZA FINAL DE MEMÓRIA ---
    free_cff_matrix(cff_old_new, rows_cff1);
    free_cff_matrix(cff_new_old, rows_cff2);
    free_cff_matrix(cff_new, rows_cff3);
    free_evaluation_matrix(evals_old, num_points, num_old_polys, ctx);
    free_evaluation_matrix(evals_new, num_points, num_new_polys, ctx);
    for (long i = 0; i < combos.capacity_old; i++) { fq_nmod_clear(combos.combos_old[i].x, ctx); fq_nmod_clear(combos.combos_old[i].y, ctx); }
    free(combos.combos_old);
    for (long i = 0; i < combos.capacity_new; i++) { fq_nmod_clear(combos.combos_new[i].x, ctx); fq_nmod_clear(combos.combos_new[i].y, ctx); }
    free(combos.combos_new);
    if (old_polys != NULL) { for (long i = 0; i < num_old_polys; i++) fq_nmod_poly_clear(old_polys[i], ctx); free(old_polys); }
    if (new_polys != NULL) { for (long i = 0; i < num_new_polys; i++) fq_nmod_poly_clear(new_polys[i], ctx); free(new_polys); }
    if (all_polys != NULL) { for (long i = 0; i < num_all_polys; i++) fq_nmod_poly_clear(all_polys[i], ctx); free(all_polys); }
    if (partitions != NULL) { for (int i = 0; i < num_steps; i++) { for (long j = 0; j < partitions[i].capacity_all; j++) fq_nmod_clear(partitions[i].all_elements[j], ctx); for (long j = 0; j < partitions[i].capacity_only; j++) fq_nmod_clear(partitions[i].only_elements[j], ctx); free(partitions[i].all_elements); free(partitions[i].only_elements); } free(partitions); }
    fmpz_clear(pz);
    fq_nmod_ctx_clear(ctx);
    
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

        // AQUI ESTÁ A LÓGICA PRINCIPAL:
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

/**
 * @brief Libera a memória de uma matriz CFF (int**) alocada dinamicamente.
 */
void free_cff_matrix(int** matrix, long num_rows) {
    if (matrix == NULL) return;

    for (long i = 0; i < num_rows; i++) {
        free(matrix[i]); // Libera cada linha
    }
    free(matrix); // Libera o ponteiro principal
}

/**
 * @brief Imprime uma matriz CFF (int**).
 * Limita a impressão a um canto de 8x8 para manter a saída legível.
 */
void print_cff_matrix(int** matrix, long num_rows, long num_cols) {

    for (long i = 0; i < num_rows; i++) {
        for (long j = 0; j < num_cols; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}