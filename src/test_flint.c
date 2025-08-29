// Run command -> gcc -I/opt/homebrew/include -L/opt/homebrew/lib test_flint.c -o test_flint -lflint -lgmp
#include <stdio.h>
#include <stdlib.h>
#include <string.h> // Para calloc
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fq_nmod.h"

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


int main() {
    printf("FLINT version: %s\n\n", FLINT_VERSION);

    // --- Parâmetros ---
    fmpz_t pz;
    fmpz_init(pz);
    fmpz_set_ui(pz, 2);
    long D = 4; // Campo principal será F_q^d
    
    // Subcorpos que queremos encontrar dentro de F_q
    long Fq_steps[] = {2,4}; 
    int num_steps = 2;

    // --- Inicialização ---
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, pz, D, "a");
    printf("Contexto inicializado para F_%ld^%ld.\n", fmpz_get_ui(pz), D);
    printf("Particionando pelos subcorpos de tamanho: %ld\n\n", Fq_steps[0]);
    
    // --- Execução da Lógica de Particionamento ---
    subfield_partition* results = partition_by_subfields(Fq_steps, num_steps, ctx);

    // --- Impressão dos Resultados ---
    for (int i = 0; i < num_steps; i++) {
        printf("--- Subcorpo F_%ld ---\n", results[i].q);
        printf("  Total de elementos (F_sets): %ld\n", results[i].count_all);
        printf("  Elementos 'Only': %ld\n", results[i].count_only);
        
        printf("  Conteúdo de 'All': { ");
        for(long j = 0; j < results[i].count_all; j++) {
            fq_nmod_print_pretty(results[i].all_elements[j], ctx);
            if (j < results[i].count_all - 1) printf(", ");
        }
        printf(" }\n");

        printf("  Conteúdo de 'Only': { ");
        for(long j = 0; j < results[i].count_only; j++) {
            fq_nmod_print_pretty(results[i].only_elements[j], ctx);
            if (j < results[i].count_only - 1) printf(", ");
        }
        printf(" }\n\n");
    }

    // --- Liberação de Memória (MUITO IMPORTANTE) ---
    printf("Limpando a memória...\n");
    for (int i = 0; i < num_steps; i++) {
        for (long j = 0; j < results[i].capacity_all; j++) fq_nmod_clear(results[i].all_elements[j], ctx);
        for (long j = 0; j < results[i].capacity_only; j++) fq_nmod_clear(results[i].only_elements[j], ctx);
        free(results[i].all_elements);
        free(results[i].only_elements);
    }
    free(results);
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
 * @brief Workaround para o bug em fq_nmod_set_ui. Constrói o elemento via aritmética.
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