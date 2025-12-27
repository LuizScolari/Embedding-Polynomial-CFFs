/**
 * @file main_benchmark.c
 * @brief Entry point for Cover-Free Families (CFFs) generation with benchmark support.
 * 
 * This program allows generating initial CFFs or expanding existing CFFs
 * using polynomial or monotone constructions over finite fields.
 * 
 * BENCHMARK VERSION: Includes automated timing measurements.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cff_builder_benchmark.h"
#include <sys/stat.h>

/* 
 *  BENCHMARK CONFIGURATION
 *  
 *  Change BENCHMARK_ITERATIONS to reduce execution time.
 *  Default: 100 iterations per test.
 */

/** @brief Number of iterations for benchmark averaging */
#define BENCHMARK_ITERATIONS 100

/* External benchmark variables from cff_builder.c */
extern int benchmark_mode;
extern int benchmark_last_iteration;
extern double benchmark_time1_accumulated;
extern double benchmark_time2_accumulated;

/**
 * @brief Resets all benchmark accumulators.
 */
void reset_benchmark_accumulators(void) {
    benchmark_time1_accumulated = 0.0;
    benchmark_time2_accumulated = 0.0;
}

/**
 * @brief Runs a single benchmark test for 'f' action (generate from scratch).
 * 
 * @param construction Construction type ('p' or 'm').
 * @param fq Finite field size.
 * @param k Maximum polynomial degree.
 */
void run_benchmark_f(char construction, long fq, long k) {
    printf("\n");
    printf("================================================================================\n");
    printf("  TEST: ./generate_cff %c f %ld %ld\n", construction, fq, k);
    printf("  Running %d iterations...\n", BENCHMARK_ITERATIONS);
    printf("================================================================================\n");
    
    reset_benchmark_accumulators();
    benchmark_mode = 1;
    
    for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
        /* Mark last iteration for file saving */
        benchmark_last_iteration = (i == BENCHMARK_ITERATIONS - 1) ? 1 : 0;
        
        generate_cff(construction, fq, k);
    }
    
    benchmark_mode = 0;
    benchmark_last_iteration = 0;
    
    /* Calculate averages */
    double avg_time1 = benchmark_time1_accumulated / BENCHMARK_ITERATIONS;
    double avg_time2 = benchmark_time2_accumulated / BENCHMARK_ITERATIONS;
    
    printf("\n");
    printf("  RESULTS:\n");
    printf("    Iterations: %d\n", BENCHMARK_ITERATIONS);
    printf("    Time 1 (Inverted Index + CFF Generation):  %.6f seconds (average)\n", avg_time1);
    printf("    Time 2 (Only CFF Matrix Generation):       %.6f seconds (average)\n", avg_time2);
    printf("================================================================================\n");
}

/**
 * @brief Runs a single benchmark test for 'g' action (embedding).
 * 
 * @param construction Construction type ('p' or 'm').
 * @param d CFF parameter d (for monotone).
 * @param Fq_steps Array with finite field sizes.
 * @param k_steps Array with maximum polynomial degrees.
 */
void run_benchmark_g(char construction, int d, long* Fq_steps, long* k_steps) {
    printf("\n");
    printf("================================================================================\n");
    if (construction == 'p') {
        printf("  TEST: ./generate_cff_benchmark %c g %ld %ld %ld %ld\n", 
               construction, Fq_steps[0], Fq_steps[1], k_steps[0], k_steps[1]);
    } else {
        printf("  TEST: ./generate_cff_benchmark %c g %d %ld %ld %ld %ld\n", 
               construction, d, Fq_steps[0], Fq_steps[1], k_steps[0], k_steps[1]);
    }
    printf("  Running %d iterations...\n", BENCHMARK_ITERATIONS);
    printf("================================================================================\n");
    
    reset_benchmark_accumulators();
    benchmark_mode = 1;
    
    for (int i = 0; i < BENCHMARK_ITERATIONS; i++) {
        /* Mark last iteration for file saving */
        benchmark_last_iteration = (i == BENCHMARK_ITERATIONS - 1) ? 1 : 0;
        
        embed_cff(construction, d, Fq_steps, k_steps);
    }
    
    benchmark_mode = 0;
    benchmark_last_iteration = 0;
    
    /* Calculate averages */
    double avg_time1 = benchmark_time1_accumulated / BENCHMARK_ITERATIONS;
    double avg_time2 = benchmark_time2_accumulated / BENCHMARK_ITERATIONS;
    
    printf("\n");
    printf("  RESULTS:\n");
    printf("    Iterations: %d\n", BENCHMARK_ITERATIONS);
    printf("    Time 1 (Inverted Index + CFF Generation + Concatenation): %.6f seconds (average)\n", avg_time1);
    printf("    Time 2 (Only CFF Matrix Generation + Concatenation):      %.6f seconds (average)\n", avg_time2);
    printf("================================================================================\n");
}

/**
 * @brief Prints usage information.
 */
void print_usage(const char* program_name) {
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  Normal mode:\n");
    fprintf(stderr, "    Initial CFF:       %s p f <q> <k>\n", program_name);
    fprintf(stderr, "    Embedding CFFs:    %s p g <q0> <q1> <k0> <k1>\n", program_name);
    fprintf(stderr, "    Monotone CFFs:     %s m g <d> <q0> <q1> <k0> <k1>\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "  Benchmark mode:\n");
    fprintf(stderr, "    Single benchmark 'f':    %s benchmark p f <q> <k>\n", program_name);
    fprintf(stderr, "    Single benchmark 'g':    %s benchmark p g <q0> <q1> <k0> <k1>\n", program_name);
    fprintf(stderr, "    Single benchmark 'g' m:  %s benchmark m g <d> <q0> <q1> <k0> <k1>\n", program_name);
    fprintf(stderr, "\n");
    fprintf(stderr, "  Use run_benchmarks.sh for automated test suites (F2, F3, F5, MONOTONE).\n");
}

/**
 * @brief Main function of the program.
 * 
 * Processes command line arguments to generate or expand CFFs,
 * with optional benchmark mode.
 * 
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @return 0 on success, 1 on error.
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    mkdir("CFFs", 0777);

    /* Check for benchmark mode */
    if (strcmp(argv[1], "benchmark") == 0) {
        /* Single benchmark with specific parameters */
        if (argc < 4) {
            fprintf(stderr, "Error: Insufficient arguments for benchmark.\n");
            print_usage(argv[0]);
            return 1;
        }
        
        char construction = argv[2][0];
        char action = argv[3][0];
        
        if (action == 'f') {
            if (construction != 'p') {
                fprintf(stderr, "Error: Action 'f' is only available for polynomial construction ('p').\n");
                return 1;
            }
            if (argc != 6) {
                fprintf(stderr, "Error: Benchmark 'p f' requires: benchmark p f <q> <k>\n");
                return 1;
            }
            
            long fq = atol(argv[4]);
            long k = atol(argv[5]);
            run_benchmark_f(construction, fq, k);
            
        } else if (action == 'g') {
            long* Fq_steps = malloc(2 * sizeof(long));
            long* k_steps = malloc(2 * sizeof(long));
            int d = 0;
            
            if (construction == 'p') {
                if (argc != 8) {
                    fprintf(stderr, "Error: Benchmark 'p g' requires: benchmark p g <q0> <q1> <k0> <k1>\n");
                    free(Fq_steps); free(k_steps);
                    return 1;
                }

                Fq_steps[0] = atol(argv[4]); 
                Fq_steps[1] = atol(argv[5]);
                k_steps[0] = atol(argv[6]);  
                k_steps[1] = atol(argv[7]);
            } else if (construction == 'm') {
                if (argc != 9) {
                    fprintf(stderr, "Error: Benchmark 'm g' requires: benchmark m g <d> <q0> <q1> <k0> <k1>\n");
                    free(Fq_steps); free(k_steps);
                    return 1;
                }
                d = atol(argv[4]); 
                Fq_steps[0] = atol(argv[5]); 
                Fq_steps[1] = atol(argv[6]);
                k_steps[0] = atol(argv[7]);  
                k_steps[1] = atol(argv[8]);
            } else {
                fprintf(stderr, "Error: Unknown construction type '%c'.\n", construction);
                free(Fq_steps); free(k_steps);
                return 1;
            }
            
            run_benchmark_g(construction, d, Fq_steps, k_steps);
            
            free(Fq_steps);
            free(k_steps);
        } else {
            fprintf(stderr, "Error: Unknown action '%c' for benchmark. Use 'f' or 'g'.\n", action);
            return 1;
        }
        
        return 0;
    }

    /* Normal mode (non-benchmark) */
    if (argc < 3) {
        fprintf(stderr, "Error: Insufficient arguments.\n");
        print_usage(argv[0]);
        return 1;
    }

    char construction = argv[1][0];
    char action = argv[2][0];

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

        embed_cff(construction, d, Fq_steps, k_steps);

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

        generate_cff(construction, fq, k);

    } else {
        fprintf(stderr, "Error: Unknown action '%c'. Use 'g' or 'f'.\n", action);
        return 1;
    }

    return 0;
}