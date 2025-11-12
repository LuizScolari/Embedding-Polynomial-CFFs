#ifndef CFF_FILE_GENERATOR_H
#define CFF_FILE_GENERATOR_H

struct cff_parameters {
    char construction;
    int* Fqs;
    int fqs_count;
    int* ks;
    int ks_count; 
};

int** read_cff_from_file(const char* filename, long* rows, long* cols);
void write_cff_to_file(const char* filename, char construction, long* Fq_steps, int fqs_count, long* K_steps, int ks_count,  int** matrix, long rows, long cols);
struct cff_parameters* read_parameters(const char* filename);

#endif 