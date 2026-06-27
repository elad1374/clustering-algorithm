# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include <stdlib.h>
# include "symnmf.h"


/* convert 2 dim Python lists to 2 dim double arrays */
/* returns new converted array or null on failure */
/* Py_list - address of python 2 dim list, rows - amount of lists in the 2 dim list, cols - amount of values in every list */
double** python_lists_to_arrays(PyObject* Py_list, int rows, int cols) {
    
    PyObject *curr_row;
    double *new_row;
    double **result;
    int i, j;

    result = (double**)malloc(sizeof(double*) * rows); /* create head array */

    for(i = 0; i < rows; i++){ /* create arrays */
        curr_row = PyList_GetItem(Py_list, i); /* get row */
        if(curr_row == NULL) {
            free_matrix(result, i -1);
            return NULL;
        }

        new_row = (double*)malloc(sizeof(double) * cols); /* create array */
        if(new_row == NULL) {
            free_matrix(result, i -1);
            return NULL;
        }

        for(j = 0; j < cols; j++) /* attach array's values */
            new_row[j] = PyFloat_AsDouble(PyList_GetItem(curr_row, j));

        result[i] = new_row; /* attach arrays */
    }

    return result;
}

/* convert double 2 dim arrays to 2 dim Python lists */
/* returns new converted list or null on failure */
/* matrix - address of 2 dim array, rows - amount of arrays in the 2 dim array, cols - amount of values in every array */
PyObject *arrays_to_python_lists(double **matrix, int rows, int cols){
    
    PyObject *new_row, *result, *python_double;
    int i, j;

    if(!(result = PyList_New(rows))) /* create head list */
        return NULL;

    for(i = 0; i < rows; i++){ /* create lists */
        if(!(new_row = PyList_New(cols)))
            return NULL;

        for(j = 0; j < cols; j++){ /* attach list's values */
            if(!(python_double = Py_BuildValue("d",matrix[i][j]))) 
                return NULL;
            
            if(PyList_SetItem(new_row, j, python_double) < 0) /* attach value */
                return NULL;
        }

        if(PyList_SetItem(result, i, new_row) < 0) /* attach list */  
            return NULL;
    }

    return result;
}


/* wrapper for sym */
/* returns the similarity A matrix or null on failure */
static PyObject* sym(PyObject* self, PyObject* args){
    
    PyObject *vectors_list, *Py_matrix, *lst;
    int rows, cols;
    double **head, **matrix;

    /* load data */
    if (!PyArg_ParseTuple(args, "O",&vectors_list)) 
        return NULL;
    
    rows = PyObject_Length(vectors_list); /* save length of rows */
    if (rows <= 0) 
        return NULL;
    
    if (!(lst = PyList_GetItem(vectors_list, 0))) 
        return NULL;
    
    cols = PyObject_Length(lst); /* save length of columns */
    if (cols <= 0) 
        return NULL;

    if(!(head = python_lists_to_arrays(vectors_list, rows, cols))) 
        return NULL;
    
    /* call symnmf.c */
    if(!(matrix = build_similarity_matrix(head, rows, cols))){
        free_matrix(head, rows);
        return NULL;
    }
    
    /* return value */
    Py_matrix = arrays_to_python_lists(matrix, rows, rows); /* will be NULL if failed */
    
    free_matrix(matrix, rows);
    free_matrix(head, rows);
    return Py_matrix; 
}

/* wrapper for ddg */
/* returns the diagonal degree D matrix or null on failure */
static PyObject* ddg(PyObject* self, PyObject* args) {
    
    PyObject *vectors_list, *Py_matrix, *lst;
    int rows, cols;
    double **head, **matrix;
    
    /* load data */
    if (!PyArg_ParseTuple(args, "O",&vectors_list)) 
        return NULL;
    
    rows = PyObject_Length(vectors_list); /* save length of rows */
    if (rows <= 0) 
        return NULL;
    
    if (!(lst = PyList_GetItem(vectors_list, 0))) 
        return NULL;
    
    cols = PyObject_Length(lst); /* save length of columns */
    if (cols <= 0) 
        return NULL;

    if(!(head = python_lists_to_arrays(vectors_list, rows, cols))) 
        return NULL;
    
    /* call symnmf.c */
    if(!(matrix = ddg_func(head, rows, cols))){
        free_matrix(head, rows);
        return NULL;
    }
    
    /* return value */
    Py_matrix = arrays_to_python_lists(matrix, rows, rows); /* will be NULL if failed */
    
    free_matrix(matrix, rows);
    free_matrix(head, rows);
    return Py_matrix;
}

/* wrapper for norm */
/* returns the normalized similarity W matrix or null on failure */
static PyObject* norm(PyObject* self, PyObject* args) {
    
    PyObject *vectors_list, *Py_matrix, *lst;
    int rows, cols;
    double **head, **matrix;
    
    /* load data */
    if (!PyArg_ParseTuple(args, "O",&vectors_list)) 
        return NULL;
    
    rows = PyObject_Length(vectors_list); /* save length of rows */
    if (rows <= 0) 
        return NULL;
    
    if (!(lst = PyList_GetItem(vectors_list, 0))) 
        return NULL;
    
    cols = PyObject_Length(lst); /* save length of columns */
    if (cols <= 0) 
        return NULL;

    if(!(head = python_lists_to_arrays(vectors_list, rows, cols))) 
        return NULL;
    
    /* call symnmf.c */
    if(!(matrix = norm_func(head, rows, cols))){
        free_matrix(head, rows);
        return NULL;
    }
    
    /* return value */
    Py_matrix = arrays_to_python_lists(matrix, rows, rows); /* will be NULL if failed */
    
    free_matrix(matrix, rows);
    free_matrix(head, rows);
    return Py_matrix; 
}

/* wrapper for symnmf */
/* returns the last H matrix or null on failure */
static PyObject* symnmf(PyObject* self, PyObject* args) {
    PyObject *Py_initial_matrix_H, *Py_matrix_H, *Py_matrix_W, *lst;
    int N_size, K;
    double **initial_matrix_H, **matrix_W, **final_matrxi_H;
    /* load data */
    if (!PyArg_ParseTuple(args,"OO",&Py_initial_matrix_H,&Py_matrix_W)) 
        return NULL;
    
    N_size = PyObject_Length(Py_initial_matrix_H); /* save size of matrix */
    if (N_size <= 0) 
        return NULL;
    if (!(lst = PyList_GetItem(Py_initial_matrix_H, 0))) 
        return NULL;

    K = PyObject_Length(lst); /* save K size*/
    if (K <= 0) 
        return NULL;
    if(!(initial_matrix_H = python_lists_to_arrays(Py_initial_matrix_H, N_size, K))) 
        return NULL;
    if(!(matrix_W = python_lists_to_arrays(Py_matrix_W, N_size, N_size))){
        free_matrix(initial_matrix_H, N_size);
        return NULL;
    }

    /* call symnmf.c */
    final_matrxi_H = symnmf_func(initial_matrix_H, matrix_W, N_size, K);
    free_matrix(initial_matrix_H, N_size);
    free_matrix(matrix_W, N_size);

    /* return result */
    if(final_matrxi_H == NULL)
        return NULL;

    Py_matrix_H = arrays_to_python_lists(final_matrxi_H, N_size, K); /* will be NULL if failed */
    free_matrix(final_matrxi_H, N_size);
    return Py_matrix_H; 
}


/* API with Python */

static PyMethodDef symnmfMethods[] = {
    {
        "sym",                   
        (PyCFunction) sym, 
        METH_VARARGS,           
        PyDoc_STR("returns the similarity A matrix")    
    }, {
        "ddg",
        (PyCFunction) ddg,
        METH_VARARGS,
        PyDoc_STR("returns the diagonal degree D matrix")
    }, {
        "norm",
        (PyCFunction) norm,
        METH_VARARGS,
        PyDoc_STR("returns the normalized similarity W matrix")
    }, {
        "symnmf",
        (PyCFunction) symnmf,
        METH_VARARGS,
        PyDoc_STR("returns the last H matrix")
    }, 
    {NULL, NULL, 0, NULL}     
};

static struct PyModuleDef symnmf_module = {
    PyModuleDef_HEAD_INIT,
    "symnmf", 
    NULL, 
    -1,  
    symnmfMethods 
};

PyMODINIT_FUNC PyInit_symnmf(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmf_module);
    if (!m) {
        return NULL;
    }
    return m;
}