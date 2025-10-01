#include "cff_builder.h"

int main() {
    // Parâmetros para a nova extensão (ex: ir de F_4 para F_16)
    long Fq_steps[] = {2,4};
    long k_steps[] = {1,2};
    int num_steps = 2;
    
    embeed_cff(Fq_steps, k_steps, num_steps);

    return 0;
}