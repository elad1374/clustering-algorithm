# include <stdio.h> 
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include "prototypes.h"


int N = -1; /* size of matrix N*N */
int K = -1; /* size of K centroids */

/* constants */
char *ERROR = "An Error Has Occurred\n";
int MAX_ITER = 300;
double EPSILON = 1e-4;
double ADDER = 1e-6;
double BETA = 0.5;

/* main program */
/* gets desired goal and a file path to file with a list of vectors and prints the desired matrix's goal */
/* first value is the goal and second value is the file path */
/* prints "An Error Has Occurred" on any error */
int main(int argc, char *argv[]){
    char *goal, *file_name;
    double **matrix, **vector_list;
    int rows, cols;

    /* read inputs */
    if(argc != 3) 
        failed(); 
    goal = argv[1];
    file_name = argv[2];
    
    if(strcmp(goal, "sym") != 0 && strcmp(goal, "ddg") != 0 && strcmp(goal, "norm") != 0) /* validate goal */
        failed();
    get_dimentions(file_name, &rows, &cols);
    if(rows <= 0 || cols <= 0)  /* validate rows and cols */
        failed();

    if(!(vector_list = get_vectors_from_file(file_name, rows, cols))) /* get vector list from file */
        failed();

    /* classify the goal and carry it out */
    if(strcmp(goal, "sym") == 0) 
        matrix = build_similarity_matrix(vector_list, rows, cols);
    else if(strcmp(goal, "ddg") == 0) 
        matrix = ddg_func(vector_list, rows, cols);
    else /* goal == norm */
        matrix = norm_func(vector_list, rows, cols);
    free_matrix(vector_list, rows);

    if(matrix == NULL) /* failure, print error massage and exit */
        failed();

    print_matrix(matrix, rows, rows); /* print result */
    free_matrix(matrix, rows);
    return 0;
}

/* program failed - prints ERROR massage and exit the program */
/* pre - all memory allocations are freed */
void failed(){
    printf("%s",ERROR);
    exit(1);
}

/* returns the similarity A matrix or null on failure */
/* list - address to list of vectors(2 dim double array), rows - amount of vectors, cols - amount of vector's values */
double **build_similarity_matrix(double **list, int rows, int cols){
    
    int i, j;
    double **matrix_A, *new_row, value;
    N = rows; /* update size of matrix */

    if(!(matrix_A = (double**)malloc(sizeof(double*) * N))) /* create head array */
        return NULL;
    
    for(i = 0; i < N; i++)
    {   
        new_row = (double*)malloc(sizeof(double) * N);
        if(new_row == NULL){
            free_matrix(matrix_A, i -1);
            return NULL;
        }

        for(j = 0; j < N; j++)
        {
            value = 0; /* calculate value */
            if(i != j) 
                value = similarity_matrix_math(list[i], list[j], cols);
            
            if(value < 0){ /* check for seccesful value calcuation */
                free(new_row);
                free_matrix(matrix_A, i -1);
                return NULL;
            }
            new_row[j] = value;
        }

        matrix_A[i] = new_row; /* attach row */
    }

    return matrix_A;
}

/* returns the diagonal degree D matrix or null on failure */
/* list - address to list of vectors(2 dim double array), rows - amount of vectors, cols - amount of vector's values */
double **ddg_func(double **list, int rows, int cols){
    
    double **matrix_A, **matrix_D;
    N = rows; /* update size of matrix */

    if(!(matrix_A = build_similarity_matrix(list, rows, cols))) 
        return NULL;
    matrix_D = build_diagonal_degree_matrix(matrix_A);

    free_matrix(matrix_A, N);
    return matrix_D;
}

/* returns the normalized similarity W matrix or null on failure */
/* list - address to list of vectors(2 dim double array), rows - amount of vectors, cols - amount of vector's values */
double** norm_func(double **list, int rows, int cols){

    double **matrix_A, **matrix_D, **matrix_W;
    N = rows; /* update size of matrix */

    if(!(matrix_A = build_similarity_matrix(list, rows, cols))) 
        return NULL;
    if(!(matrix_D = build_diagonal_degree_matrix(matrix_A))){
        free_matrix(matrix_A, N);
        return NULL;
    }
    
    matrix_W = build_normalized_similarity_matrix(matrix_A, matrix_D);
    
    free_matrix(matrix_A, N);
    free_matrix(matrix_D, N);
    return matrix_W;
}

/* returns the last H matrix or null on failure */
/* initial_H_matrix - address of the first random computed H matrix, matrix_W - address of the computed normalized similarity matrix */
/* N_size/K_size - dimentions of H matrix N*K */
double **symnmf_func(double **initial_H_matrix, double **matrix_W, int N_size, int K_size){
    
    double **new_H_matrix = NULL, **old_H_matrix;
    int iterations_left = MAX_ITER;
    N = N_size; /* update size of matrix */
    K = K_size;

    old_H_matrix = initial_H_matrix;
    while(iterations_left--)
    {
        new_H_matrix = get_updated_H(old_H_matrix, matrix_W);
        
        if(new_H_matrix == NULL){
            if(old_H_matrix != initial_H_matrix)
                free_matrix(old_H_matrix, N);
            return NULL;
        }
        
        if(euclidean_distances_squared_matrix(new_H_matrix ,old_H_matrix, N, K) < EPSILON){
            if(old_H_matrix != initial_H_matrix)
                free_matrix(old_H_matrix, N);
            return new_H_matrix;
        }
        
        if(old_H_matrix != initial_H_matrix)
            free_matrix(old_H_matrix, N);
        old_H_matrix = new_H_matrix;
    }
    
    return new_H_matrix;
}

/* returns matrix H^(t) based on matrix H^(t-1) and matrix W or null on failure */
/* matrix_H - address of the current computed H matrix, matrix_W - address of the computed normalized similarity matrix */
double **get_updated_H(double **matrix_H, double **matrix_W){
    int i, j;
    double **new_matrix, **upper_calc, **lower_calc;

    if(!(upper_calc = mul_matrix(matrix_W, matrix_H, N, K, N))) /* calculate upper part for H */
        return NULL;
    if(!(lower_calc = calc_lower_part_of_updated_H(matrix_H))){ 
        free_matrix(upper_calc, N);
        return NULL;
    }

    /* build the new H */    
    if(!(new_matrix = (double**)malloc(sizeof(double*) * N))){
        free_matrix(upper_calc, N);
        free_matrix(lower_calc, N);
        return NULL;
    }
    
    for(i = 0; i < N; i++){
        if(!(new_matrix[i] = (double*)malloc(sizeof(double) * K))){ /* allocate space for new row */
            free_matrix(new_matrix, i -1);
            free_matrix(upper_calc, N);
            free_matrix(lower_calc, N);
            return NULL;
        }
        
        for(j = 0; j < K; j++){ /* apply row's value*/
            /* calculate value */
            if(lower_calc[i][j] != 0) 
                new_matrix[i][j] = matrix_H[i][j] * (1 - BETA + BETA * upper_calc[i][j]/lower_calc[i][j]);
            else /* won't divide by 0 */
                new_matrix[i][j] = matrix_H[i][j] * (1 - BETA + BETA * upper_calc[i][j]/ADDER);
        }
    }
    free_matrix(upper_calc, N);
    free_matrix(lower_calc, N);
    return new_matrix;
}

/* returns the lower part needed for calculating the new matrix H */
/* matrix_H - address of the old H matrix */
double **calc_lower_part_of_updated_H(double **matrix_H){

    double **lower_calc, **tmp;

    if(!(tmp = matrix_transposed(matrix_H, N, K))) /* transpose matrix H */
        return NULL;
    
    lower_calc = mul_matrix(matrix_H, tmp, N, N, K);/* initial lower calc - first 2 matrix */
    free_matrix(tmp, K);
    if(lower_calc == NULL)
        return NULL;

    tmp = lower_calc;
    lower_calc = mul_matrix(lower_calc, matrix_H, N, K, N);/* final calc - last matrix */
    free_matrix(tmp, N);
    if(lower_calc == NULL)
        return NULL;
    
    return lower_calc;
}

/* returns the diagonal degree matrix D based on the similarity matrix or null on failure */
/* matrix_A - adress of the similarity matrix A */
double **build_diagonal_degree_matrix(double **matrix_A){

    int i;
    double **matrix_D;

    if(!(matrix_D = (double**)malloc(sizeof(double*) * N))) /* allocate memory for D matrix */
        return NULL;
        
    for(i = 0; i < N; i++)
    {
        if(!(matrix_D[i] = (double*)calloc(N, sizeof(double)))){ /* create rows */
            free_matrix(matrix_D, i -1);
            return NULL;
        }

        matrix_D[i][i] = sum_row(matrix_A[i], N); /* apply value */
    }

    return matrix_D;
}

/* returns the normalized similarity matrix W based on the similarity matrix A */
/* and the digonal degree matrix D or null on failure */
/* matrix_A - adress of the similarity matrix A, matrix_D - adress of the digonal degree matrix D */
double** build_normalized_similarity_matrix(double **matrix_A, double **matrix_D){

    double **matrix_W, **first_calc;
    int i;
    
    if(matrix_A == NULL || matrix_D == NULL) 
        return NULL;
    
    for(i = 0; i < N; i++) /* modify matrix D */
        matrix_D[i][i] = 1 / sqrt(matrix_D[i][i]);
    
    /* calculate W */    
    if(!(first_calc = mul_matrix(matrix_D, matrix_A, N, N, N)))
        return NULL;
    matrix_W = mul_matrix(first_calc, matrix_D, N, N, N);

    free_matrix(first_calc, N);
    return matrix_W;
}

/* returns the calculated value of an element in the similarity matrix A of the 2 vectors or null on failure */
/* first - adress of the first vector(double array), second - adress of the second vector, cols - length of the arrays */
double similarity_matrix_math(double *first, double *second, int cols){
    if(first == NULL || second == NULL) return -1;
    return exp(-(euclidean_distances_squared_vector(first, second, cols)/2));
}

/* returns the euclidean distances squared vector(double array) of 2 vectors or -1 on failure */
/* first - adress of the first vector(double array), second - adress of the second vector, cols - length of the arrays */
double euclidean_distances_squared_vector(double *first, double *second, int cols){
    
    double sum = 0.0;
    int i;

    if(first == NULL || second == NULL) 
        return -1;

    for(i = 0; i < cols; i++)
        sum += pow(first[i] - second[i], 2);
    
    return sum;
}

/* returns the multiplicity of 2 matrix(2 dim double arrays) or null on failure */
/* pre - cols_first == rows_second*/
/* first - adress of the first matrix to multiply, second - adress of the second matrix to multiply */
/* first_rows - length of rows in first matrix, second_cols - length of cols in the second matrix */
/* first_cols_and_second_rows - the length if cols/rows of the first/second matrix */
double** mul_matrix(double **first, double **second, int first_rows, int second_cols, int first_cols_and_second_rows){
    
    int i, j, r;
    double **result;

    /* create empty matrix */
    if (!(result = (double**)malloc(sizeof(double*) * first_rows))) 
        return NULL;

    for (i = 0; i < first_rows; i++){
        if (!(result[i] = (double*)malloc(sizeof(double) * second_cols))){
            free_matrix(result, i - 1);
            return NULL;
        }
    }

    /* multiply both matrices */
    for (i = 0; i < first_rows; i++) {
        for (j = 0; j < second_cols; j++) {
            result[i][j] = 0.0;
            for (r = 0; r < first_cols_and_second_rows; r++) {
                result[i][j] += first[i][r] * second[r][j];
            }
        }
    }

    return result;
}

/* returns the euclidean distances squared of matrices(2 dim double array) between 2 matrices or -1 on failure */
/* first - adress of the first matrix(2 dim double array), second - adress of the second matrix */ 
/* rows - amount of arrays in matrix, cols - length of the arrays */
double euclidean_distances_squared_matrix(double **first, double **second, int rows, int cols){
    
    double sum = 0.0;
    int i;

    if(first == NULL || second == NULL) /* check viability of matrices */
        return -1;

    for(i = 0; i < rows; i++) /* calculate */
        sum += euclidean_distances_squared_vector(first[i], second[i], cols);

    return sum;
}

/* prints a vector(double array) */
/* row - the address of the vector to print, length - the length of the vector */
void print_row(double *row, int length){
    
    int i;

    if(row == NULL || length == 0) 
        return;
    
    for(i = 0; i < length -1; i++) /* print all values except the last */
        printf("%.4f,", row[i]);
    printf("%.4f\n", row[length -1]); /* print the last value */
}

/* prints a matrix(2 dim double array) */
/* matrix - the address of the matrix, rows - amount of vectors in matrix, cols - length of the vectors in matrix */
void print_matrix(double **matrix, int rows, int cols){
    
    int i;

    if(matrix == NULL || rows == 0 || cols == 0) 
        return;
    
    for(i = 0; i < rows; i++)
        print_row(matrix[i], cols);
}

/* return a the sum of vector(double array) */
/* row - address of the vector, length - the length of the vector */
double sum_row(double *row, int length){
    double sum = 0;
    int i;
    for(i = 0; i < length; i++) 
        sum += row[i];
    return sum;
}

/* returns a transposed verison of the matrix */
/* matrix - the address of the matrix, rows - amount of vectors in matrix, cols - length of the vectors in matrix */
double **matrix_transposed(double **matrix, int rows, int cols){
    
    double **transposed;
    int i, j;

    if(!(transposed = (double**)malloc(sizeof(double*) * cols))) 
        return NULL;

    for(i = 0; i < cols; i++)
    {
        if(!(transposed[i] = (double*)malloc(sizeof(double) * rows))){ /* allocate space */
            free_matrix(transposed, i - 1);
            return NULL;
        }

        for(j = 0; j < rows; j++) /* copy values */
            transposed[i][j] = matrix[j][i];
    }

    return transposed;
}

/* apply dimentions of the values in file or -1 on failure */
/* file_name - file path of the file , rows/cols - address to save the amount lines/values of the file */
void get_dimentions(char *file_name, int *rows, int *cols){
    
    FILE *file_pointer = NULL;
    char c;
    double tmp;

    /* first initialize */
    *rows = -1;
    *cols = -1;

    if(file_name == NULL)
        return;
    
    if(!(file_pointer = fopen(file_name, "r")))
        return;
    
    /* start counting */
    *rows = 0;
    *cols = 0;

    while(fscanf(file_pointer, "%lf%c", &tmp, &c) == 2)
    {
        if(c == EOF)
            break;
        
        if((*rows) == 0) /* count columns(of the first row) */
            (*cols)++;

        if(c == '\n') /* count rows */
            (*rows)++;
    }
    
    fclose(file_pointer);
}

/* gather the vectors from file and returns in as 2 dim double array or null on failure */
/* file_name - file path of the file, rows - amount of vectors in the file, cols - amount of values in every line */
double **get_vectors_from_file(char *file_name, int rows, int cols){

    FILE *file_pointer;
    double **vectors;
    char c;
    int i, j;

    /* create lists */
    if(!(vectors = (double**)malloc(sizeof(double*) * rows)))
        return NULL;
    
    /* open file */
    file_pointer = fopen(file_name, "r");
    if (file_pointer == NULL){
        free(vectors);
        return NULL;
    }

    /* read the data */
    for(i = 0; i < rows; i++)
    {
        if(!(vectors[i] = (double*)malloc(sizeof(double) * cols))){ /* allocate new row memory */
            free_matrix(vectors, i -1);
            fclose(file_pointer);
            return NULL;
        }
        for(j = 0; j < cols; j++) /* read line */
            if(fscanf(file_pointer, "%lf%c", &vectors[i][j], &c) != 2){
                free_matrix(vectors, i);
                fclose(file_pointer);
                return NULL;
            }
    }
    
    fclose(file_pointer);
    return vectors;
}

/* free 2 dim double array */
/* matrix - address of 2 dim double array , rows - amount of arrays in the 2 dim array*/
void free_matrix(double** matrix, int rows){
    int i;
    if (matrix != NULL) {
        for (i = 0; i < rows; i++) {
            free(matrix[i]);
        }
        free(matrix);
    }
}

