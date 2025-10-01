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

// Função pública que a main acessa 
void embeed_cff(long* Fq_steps, long* k_steps, int num_steps);

#endif // CFF_BUILDER_H