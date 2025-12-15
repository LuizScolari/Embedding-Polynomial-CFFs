/**
 * @file cff_builder.h
 * @brief Definitions and structures for Cover-Free Families (CFFs) construction.
 * 
 * This file contains type definitions, macros, and function prototypes
 * for generating CFFs using polynomial and monotone constructions.
 */

#ifndef CFF_BUILDER_H
#define CFF_BUILDER_H

#include <stdint.h> 
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_poly.h"

/*
 * 64-BIT BITMAP MACROS
 */

/** @brief Number of bits per word (64 bits). */
#define BITS_PER_WORD 64

/** @brief Calculates the word offset for a given bit. */
#define WORD_OFFSET(b) ((b) / BITS_PER_WORD)

/** @brief Calculates the bit offset within a word. */
#define BIT_OFFSET(b)  ((b) % BITS_PER_WORD)

/** @brief Sets the bit at position col in row_ptr. */
#define SET_BIT(row_ptr, col) ((row_ptr)[WORD_OFFSET(col)] |= (1ULL << BIT_OFFSET(col)))

/** @brief Gets the bit value at position col in row_ptr. */
#define GET_BIT(row_ptr, col) (((row_ptr)[WORD_OFFSET(col)] >> BIT_OFFSET(col)) & 1ULL)

/** @brief Calculates the number of words needed to store bits. */
#define WORDS_FOR_BITS(bits) (((bits) + BITS_PER_WORD - 1) / BITS_PER_WORD)

/* 
 *  DATA STRUCTURES
 */

/**
 * @brief Structure to store a subfield element partition.
 * 
 * Contains elements belonging to a specific subfield F_q,
 * divided into "all" (all_elements) and "only new" (only_elements).
 */
typedef struct { 
    long q;                     /**< Subfield size. */
    fq_nmod_t* all_elements;    /**< All elements of the subfield. */
    long count_all;             /**< Number of elements in all_elements. */
    long capacity_all;          /**< Allocated capacity for all_elements. */
    fq_nmod_t* only_elements;   /**< Elements exclusive to this subfield. */
    long count_only;            /**< Number of elements in only_elements. */
    long capacity_only;         /**< Allocated capacity for only_elements. */
} subfield_partition;

/**
 * @brief Structure to represent an element pair (x, y).
 * 
 * Used to represent polynomial evaluation pairs.
 */
typedef struct { 
    fq_nmod_t x;    /**< First element of the pair. */
    fq_nmod_t y;    /**< Second element of the pair. */
} element_pair;

/**
 * @brief Structure to store pair combination partitions.
 * 
 * Divides pairs into "old" (from previous steps) and "new".
 */
typedef struct { 
    element_pair* combos_old;   /**< Pairs from previous steps. */
    long count_old;             /**< Number of old pairs. */
    long capacity_old;          /**< Allocated capacity for combos_old. */
    element_pair* combos_new;   /**< Pairs from current step. */
    long count_new;             /**< Number of new pairs. */
    long capacity_new;          /**< Allocated capacity for combos_new. */
} combination_partitions;

/**
 * @brief Structure to return the three new CFF blocks.
 * 
 * Contains the cff_old_new, cff_new_old, and cff_new matrices
 * generated during the embedding process.
 */
typedef struct {
    uint64_t** cff_old_new;     /**< old_new block of the CFF. */
    long rows_old_new;          /**< Number of rows in cff_old_new. */
    long cols_old_new;          /**< Number of columns in cff_old_new. */

    uint64_t** cff_new_old;     /**< new_old block of the CFF. */
    long rows_new_old;          /**< Number of rows in cff_new_old. */
    long cols_new_old;          /**< Number of columns in cff_new_old. */
    
    uint64_t** cff_new;         /**< new_new block of the CFF. */
    long rows_new;              /**< Number of rows in cff_new. */
    long cols_new;              /**< Number of columns in cff_new. */
} generated_cffs;

/**
 * @brief Structure to store polynomial partition.
 * 
 * Divides polynomials into "old" (from previous steps), 
 * "new" (exclusive to the last step), and "all".
 */
typedef struct {
    fq_nmod_poly_t* old_polys;  /**< Polynomials from previous steps. */
    long num_old_polys;         /**< Number of old polynomials. */
    fq_nmod_poly_t* new_polys;  /**< Polynomials from current step. */
    long num_new_polys;         /**< Number of new polynomials. */
    fq_nmod_poly_t* all_polys;  /**< All polynomials. */
    long num_all_polys;         /**< Total number of polynomials. */
} polynomial_partition;

/*
 * PUBLIC FUNCTION PROTOTYPES
 */

/**
 * @brief Performs embedding of an existing CFF to a larger field.
 * 
 * @param construction Construction type ('p' for polynomial, 'm' for monotone).
 * @param d CFF parameter d.
 * @param Fq_steps Array with finite field sizes.
 * @param k_steps Array with maximum polynomial degrees.
 */
void embeed_cff(char construction, int d, long* Fq_steps, long* k_steps);

/**
 * @brief Generates an initial CFF from basic parameters.
 * 
 * @param construction Construction type ('p' for polynomial, 'm' for monotone).
 * @param d CFF parameter d.
 * @param fq Finite field size.
 * @param k Maximum polynomial degree.
 */
void generate_cff(char construction, int d, long fq, long k);

#endif /* CFF_BUILDER_H */
