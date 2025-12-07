#ifndef CFF_BUILDER_H
#define CFF_BUILDER_H

#include <stdint.h> 
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_poly.h"

// --- Macros para Bitmap de 64 bits ---
#define BITS_PER_WORD 64
#define WORD_OFFSET(b) ((b) / BITS_PER_WORD)
#define BIT_OFFSET(b)  ((b) % BITS_PER_WORD)

#define SET_BIT(row_ptr, col) ((row_ptr)[WORD_OFFSET(col)] |= (1ULL << BIT_OFFSET(col)))
#define GET_BIT(row_ptr, col) (((row_ptr)[WORD_OFFSET(col)] >> BIT_OFFSET(col)) & 1ULL)

#define WORDS_FOR_BITS(bits) (((bits) + BITS_PER_WORD - 1) / BITS_PER_WORD)

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
    uint64_t** cff_old_new;
    long rows_old_new;
    long cols_old_new;

    uint64_t** cff_new_old;
    long rows_new_old;
    long cols_new_old;
    
    uint64_t** cff_new;
    long rows_new;
    long cols_new;
} generated_cffs;

// Função pública que a main acessa 
void embeed_cff(char construction, long* Fq_steps, long* k_steps);
void generate_cff(char construction, long fq, long k);

#endif // CFF_BUILDER_H