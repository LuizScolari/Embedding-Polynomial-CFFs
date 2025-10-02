#include <stdio.h>
#include <stdlib.h>
#include "cff_builder.h"
#include <sys/stat.h>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        return 1;
    }

    int num_steps = atoi(argv[1]);

    if (argc != 2 + (2 * num_steps)) {
        return 1;
    }

    long *Fq_steps = malloc(num_steps * sizeof(long));
    long *k_steps = malloc(num_steps * sizeof(long));

    if (Fq_steps == NULL || k_steps == NULL) {
        return 1;
    }

    for (int i = 0; i < num_steps; i++) {
        Fq_steps[i] = atol(argv[2 + i]);
    }

    for (int i = 0; i < num_steps; i++) {
        k_steps[i] = atol(argv[2 + num_steps + i]);
    }

    mkdir("CFFs", 0777);
    embeed_cff(Fq_steps, k_steps, num_steps);

    free(Fq_steps);
    free(k_steps);

    return 0;
}