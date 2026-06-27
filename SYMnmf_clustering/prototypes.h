
# ifndef PROTOTYPES
# define PROTOTYPES

double** build_similarity_matrix(double **list, int rows, int cols);
double** ddg_func(double **list, int rows, int cols);
double** norm_func(double **list, int rows, int cols);
double** symnmf_func(double **matrix_H, double **matrix_W, int N_size, int K_size);
void free_matrix(double** matrix, int rows);

double** build_diagonal_degree_matrix(double **matrix_A);
double** build_normalized_similarity_matrix(double **matrix_A, double **matrix_D);
double** matrix_transposed(double **matrix, int rows, int cols);
double** mul_matrix(double **first, double **second, int first_rows, int second_cols, int first_cols_and_second_rows);
double** get_updated_H(double **matrix_H, double **matrix_W);
double **calc_lower_part_of_updated_H(double **matrix_H);
double** get_vectors_from_file(char *file_name, int rows, int cols);
double similarity_matrix_math(double *first, double *second, int cols);
double euclidean_distances_squared_vector(double *first, double *second, int cols);
double euclidean_distances_squared_matrix(double **first, double **second, int rows, int cols);
double sum_row(double *row, int length);
void get_dimentions(char *file_name, int *rows, int *cols);
void print_row(double *row, int length);
void print_matrix(double **matrix, int rows, int cols);
void failed();

# endif