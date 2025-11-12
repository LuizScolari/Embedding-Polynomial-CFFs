#include <stdio.h>
#include <stdlib.h>
#include "cff_builder.h"
#include <sys/stat.h>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Erro: Argumentos insuficientes.\n");
        return 1;
    }

    char construction = argv[1][0];
    char action = argv[2][0];

    mkdir("CFFs", 0777);

    if (action == 'g') {
        if (argc != 7) {
            fprintf(stderr, "Erro (g): Número incorreto de argumentos.\n");
            return 1;
        }

        long* Fq_steps = malloc(2 * sizeof(long));
        long* k_steps = malloc(2 * sizeof(long));
        if (Fq_steps == NULL || k_steps == NULL) {
            perror("Erro: Falha ao alocar memória");
            free(Fq_steps); free(k_steps);
            return 1;
        }

        Fq_steps[0] = atol(argv[3]); 
        Fq_steps[1] = atol(argv[4]);
        k_steps[0] = atol(argv[5]);  
        k_steps[1] = atol(argv[6]);  

        embeed_cff(construction, Fq_steps, k_steps);

        free(Fq_steps);
        free(k_steps);

    } else if (action == 'f') {
        if (argc != 5) {
            fprintf(stderr, "Erro (f): Número incorreto de argumentos.\n");
            return 1;
        }

        long fq = atol(argv[3]);
        long k = atol(argv[4]);

        generate_cff(construction, fq, k);

    } else {
        fprintf(stderr, "Erro: Ação '%c' desconhecida. Use 'g' ou 'f'.\n", action);
        return 1;
    }

    return 0;
}