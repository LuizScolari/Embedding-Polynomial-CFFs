#ifndef CFF_FILE_GENERATOR_H
#define CFF_FILE_GENERATOR_H

int** read_cff_from_file(const char* filename, long* rows, long* cols);
void write_cff_to_file(const char* filename, int** matrix, long rows, long cols);

#endif 