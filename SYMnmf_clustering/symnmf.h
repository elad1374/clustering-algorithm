
# ifndef SYMNMF_H
# define SYMNMF_H

double **build_similarity_matrix(double **list, int rows, int cols);
double **ddg_func(double **list, int rows, int cols);
double **norm_func(double **list, int rows, int cols);
double **symnmf_func(double **matrix_H, double **matrix_W, int N_size, int K);

void free_matrix(double** matrix, int rows);

# endif