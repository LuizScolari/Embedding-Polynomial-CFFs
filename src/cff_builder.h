#ifndef CFF_BUILDER_H
#define CFF_BUILDER_H

#include "flint/fq_nmod.h"
#include "flint/fq_nmod_poly.h"

// --- Estruturas de Dados Públicas ---
typedef struct { 
    long q; fq_nmod_t* all_elements; 
    long count_all; 
    long capacity_all; 
    fq_nmod_t* only_elements; 
    long count_only; 
    long capacity_only; 
} subfield_partition;

typedef struct { 
    fq_nmod_t x; 
    fq_nmod_t y; 
} element_pair;

typedef struct { 
    element_pair* combos_old; 
    long count_old; 
    long capacity_old; 
    element_pair* combos_new; 
    long count_new; 
    long capacity_new; 
} combination_partitions;

// Estrutura para retornar os 3 novos blocos da CFF
typedef struct {
    int** cff_old_new;
    long rows_old_new;
    long cols_old_new;

    int** cff_new_old;
    long rows_new_old;
    long cols_new_old;
    
    int** cff_new;
    long rows_new;
    long cols_new;
} generated_cffs;


// --- Funções Públicas ---

// A função principal do nosso módulo: gera os 3 novos blocos da CFF
generated_cffs generate_new_cff_blocks(long* Fq_steps, int num_steps, long max_poly_degree);

// Funções para liberar a memória das estruturas retornadas
void free_subfield_partitions(subfield_partition* partitions, int num_steps, const fq_nmod_ctx_t ctx);
void free_combination_partitions(combination_partitions* combos, const fq_nmod_ctx_t ctx);
void free_generated_cffs(generated_cffs* cffs);

#endif // CFF_BUILDER_H