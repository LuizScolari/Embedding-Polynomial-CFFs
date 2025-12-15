/**
 * @file main.c
 * @brief Entry point for Cover-Free Families (CFFs) generation.
 * 
 * This program allows generating initial CFFs or expanding existing CFFs
 * using polynomial or monotone constructions over finite fields.
 */

#include <stdio.h>
#include <stdlib.h>
#include "cff_builder.h"
#include <sys/stat.h>

/**
 * @brief Main function of the program.
 * 
 * Processes command line arguments to generate or expand CFFs.
 * 
 * Usage:
 *   - Polynomial embedding: ./generate_cff p g <q0> <q1> <k0> <k1>
 *   - Polynomial from scratch: ./generate_cff p f <q> <k>
 *   - Monotone embedding: ./generate_cff m g <d> <q0> <q1> <k0> <k1>
 * 
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @return 0 on success, 1 on error.
 */
int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Error: Insufficient arguments.\n");
        return 1;
    }

    char construction = argv[1][0];
    char action = argv[2][0];

    mkdir("CFFs", 0777);

    if (action == 'g') {
        if (construction == 'p') {
            if (argc != 7) {
                fprintf(stderr, "Error (p g): Incorrect number of arguments.\n");
                return 1;
            }
        }

        if (construction == 'm') {
            if (argc != 8) {
                fprintf(stderr, "Error (m g): Incorrect number of arguments.\n");
                return 1;
            }
        }

        long* Fq_steps = malloc(2 * sizeof(long));
        long* k_steps = malloc(2 * sizeof(long));
        if (Fq_steps == NULL || k_steps == NULL) {
            perror("Error: Failed to allocate memory");
            free(Fq_steps); free(k_steps);
            return 1;
        }
        
        int d = 0;
        if(construction == 'p'){
            Fq_steps[0] = atol(argv[3]); 
            Fq_steps[1] = atol(argv[4]);
            k_steps[0] = atol(argv[5]);  
            k_steps[1] = atol(argv[6]); 
        } else if (construction == 'm'){
            d = atol(argv[3]); 
            Fq_steps[0] = atol(argv[4]); 
            Fq_steps[1] = atol(argv[5]);
            k_steps[0] = atol(argv[6]);  
            k_steps[1] = atol(argv[7]); 
        }

        embeed_cff(construction, d, Fq_steps, k_steps);

        free(Fq_steps);
        free(k_steps);

    } else if (action == 'f') {
        if (construction != 'p') {
            fprintf(stderr, "Error: Action 'f' is only available for polynomial construction ('p').\n");
            fprintf(stderr, "For monotone construction ('m'), use only action 'g' (embedding).\n");
            return 1;
        }
        
        if (argc != 5) {
            fprintf(stderr, "Error (p f): Incorrect number of arguments. Usage: ./generate_cff p f <q> <k>\n");
            return 1;
        }
        
        long fq = atol(argv[3]);
        long k = atol(argv[4]);

        generate_cff(construction, 0, fq, k);

    } else {
        fprintf(stderr, "Error: Unknown action '%c'. Use 'g' or 'f'.\n", action);
        return 1;
    }

    return 0;
}
