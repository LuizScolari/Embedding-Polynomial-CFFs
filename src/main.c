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
        if (construction == 'p') {
            if (argc != 7) {
                fprintf(stderr, "Erro (p g): Número incorreto de argumentos.\n");
                return 1;
            }
        }

        if (construction == 'm') {
            if (argc != 8) {
                fprintf(stderr, "Erro (m g): Número incorreto de argumentos.\n");
                return 1;
            }
        }

        long* Fq_steps = malloc(2 * sizeof(long));
        long* k_steps = malloc(2 * sizeof(long));
        if (Fq_steps == NULL || k_steps == NULL) {
            perror("Erro: Falha ao alocar memória");
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
            fprintf(stderr, "Erro: Ação 'f' só está disponível para construção polynomial ('p').\n");
            fprintf(stderr, "Para construção monotone ('m'), use apenas a ação 'g' (embedding).\n");
            return 1;
        }
        
        if (argc != 5) {
            fprintf(stderr, "Erro (p f): Número incorreto de argumentos. Uso: ./generate_cff p f <q> <k>\n");
            return 1;
        }
        
        long fq = atol(argv[3]);
        long k = atol(argv[4]);

        generate_cff(construction, 0, fq, k);

    } else {
        fprintf(stderr, "Erro: Ação '%c' desconhecida. Use 'g' ou 'f'.\n", action);
        return 1;
    }

    return 0;
}