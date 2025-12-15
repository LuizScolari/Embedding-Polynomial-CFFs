/**
 * @file cff_file_generator.h
 * @brief Definitions for CFF file reading and writing.
 * 
 * This file contains the structures and function prototypes for
 * Cover-Free Families file manipulation.
 */

#ifndef CFF_FILE_GENERATOR_H
#define CFF_FILE_GENERATOR_H

#include <stdint.h>

/*
 * DATA STRUCTURES
 */

/**
 * @brief Structure to store CFF parameters read from file.
 * 
 * Contains the construction type, d parameter, and the arrays of
 * finite field sizes and polynomial degrees.
 */
struct cff_parameters {
    char construction;  /**< Construction type ('p' or 'm'). */
    int d;              /**< CFF parameter d (for monotone construction). */
    int* Fqs;           /**< Array of finite field sizes. */
    int fqs_count;      /**< Number of elements in Fqs. */
    int* ks;            /**< Array of maximum polynomial degrees. */
    int ks_count;       /**< Number of elements in ks. */
};

/*
 * FUNCTION PROTOTYPES
 */

/**
 * @brief Reads a CFF matrix from a file.
 * 
 * @param filename Path to the file to be read.
 * @param rows Pointer to store the number of rows.
 * @param cols Pointer to store the number of columns.
 * @return CFF matrix in bitmap format.
 */
uint64_t** read_cff_from_file(const char* filename, long* rows, long* cols);

/**
 * @brief Writes a CFF matrix to a file.
 * 
 * @param filename Output file path.
 * @param construction Construction type.
 * @param d CFF parameter d.
 * @param Fq_steps Array with finite field sizes.
 * @param fqs_count Number of elements in Fq_steps.
 * @param K_steps Array with maximum polynomial degrees.
 * @param ks_count Number of elements in K_steps.
 * @param matrix CFF matrix in bitmap format.
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 */
void write_cff_to_file(const char* filename, char construction, int d, long* Fq_steps, int fqs_count, long* K_steps, int ks_count, uint64_t** matrix, long rows, long cols);

/**
 * @brief Reads CFF parameters from a file.
 * 
 * @param filename Path to the file to be read.
 * @return Pointer to structure with parameters.
 */
struct cff_parameters* read_parameters(const char* filename);

#endif /* CFF_FILE_GENERATOR_H */
