# define PY_SSIZE_T_CLEAN
# include <Python.h>

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

struct cord
{
    double value;
    struct cord *next;
};
struct vector
{
    struct vector *next;
    struct cord *cords;
};

int get_closest_cluster_index(struct vector *,struct vector *);
int cluster_not_changed(struct vector *,struct vector *);
void reset_size(int *,int );
double get_distance(struct vector *,struct vector *);
void reset_vectors(struct vector *);
void add_vectors(struct vector *,struct vector *);
void average_vectors(struct vector *,int *, struct vector *);
void free_cords(struct cord *);
void free_vector(struct vector *);
struct vector *initialize_centroids_sum(struct vector *,int);
struct vector *copy_vector_from(struct vector *);
struct vector *initialize_data(PyObject *table);


double EPSILON = 0.0;

struct vector *fit_func(struct vector *data_base,struct vector *centroids, int K, int iter, double eps)
{
    struct vector  *current;
    struct vector  *centroids_sum;
    struct vector *cur_centroid_sum;
    int *centroids_size;
    int closest; /*calculating closest centroid to vector*/
    EPSILON = eps;
    
    /*input errors*/
    if (iter <= 1 || iter >= 1000)
    {
        free_vector(data_base);
        free_vector(centroids);
        PyErr_SetString(PyExc_ValueError, "Invalid maximum iteration!\n");
        return NULL;
    }

    if(K <= 1)
    {
        free_vector(centroids);
        free_vector(data_base);
        PyErr_SetString(PyExc_ValueError, "Invalid number of clusters!\n");
        return NULL;
    }

    /*reset centroid_sum + check if database smaller then K*/
    
    centroids_sum = initialize_centroids_sum(data_base,K);
    if(centroids_sum == NULL)
    {
        free_vector(centroids);
        free_vector(data_base);
        return NULL;
    }
    
    reset_vectors(centroids_sum);

    /*array that saves the amount of vectors in centroids_sum*/
    centroids_size = calloc(K,sizeof(int));
    if(centroids_size == NULL)
    {
        free_vector(centroids_sum);
        free_vector(centroids);
        free_vector(data_base);
        return NULL;
    }
    
    while(iter--)
    {
        current = data_base;
        /*add all vectors to closest sum*/
        while (current != NULL)
        {
            
            closest = get_closest_cluster_index(current,centroids);
            
            /*add to closest sum and add +1 to its size*/
            centroids_size[closest]++;
            
            cur_centroid_sum = centroids_sum;
            
            while(closest != 0) {
                closest--;
                cur_centroid_sum = cur_centroid_sum->next;
            }
            
            add_vectors(cur_centroid_sum,current);
            
            current = current->next;
        }
        
        average_vectors(centroids_sum,centroids_size, centroids);
    
        /*swap old centroid with new one*/
        cur_centroid_sum = centroids_sum;
        centroids_sum = centroids;
        centroids = cur_centroid_sum;


        if(cluster_not_changed(centroids_sum,centroids))
            break;
        
        reset_size(centroids_size,K);
        reset_vectors(centroids_sum);
    }
    
    free_vector(centroids_sum);
    free(centroids_size);
    free_vector(data_base);
    
    return centroids;
}


/*get the index of the closest cluster*/
int get_closest_cluster_index(struct vector *vec,struct vector *centroids)
{
    
    double closest_dist = get_distance(vec,centroids),dist;
    int cnt = 0;
    int closest = 0;

    while (centroids != NULL)
    {
        dist = get_distance(vec,centroids);
        if (dist < closest_dist)
        {
            closest_dist = dist;
            closest = cnt;
        }
        centroids = centroids->next;
        cnt++;
    }
    return closest;
}

/*calculate if old and new clusters are epsilon close*/
/*new centroid need to be the same length as old centroid*/
int cluster_not_changed(struct vector *new_centroid,struct vector *old_centroid)
{
    while(new_centroid != NULL)
    {
        if(get_distance(new_centroid,old_centroid) >= EPSILON)
            return 0;

        new_centroid = new_centroid->next;
        old_centroid = old_centroid->next;
    }
    return 1;
}

/*resets an array of size k*/
void reset_size(int *centroid_size,int K)
{
    int i;
    for (i = 0;i < K;i++)
        centroid_size[i] = 0;
}

/*calculate euclidian distance between vec1 and vec2*/
double get_distance(struct vector *vec1,struct vector *vec2)
{
    double dist = 0;
    struct cord *cord1 = vec1->cords;
    struct cord *cord2 = vec2->cords;
    
    while(cord1 != NULL || cord2 != NULL)
    {
        dist += pow(cord1->value - cord2->value,2.0);
        cord1 = cord1->next;
        cord2 = cord2->next;
    }
    return pow(dist,0.5);
}

void reset_vectors(struct vector *vectors)
{
    struct cord *tmp;
    if(vectors == NULL)
        return;
    if(vectors->cords == NULL)
        return;
    

    while(vectors != NULL)
    {
        tmp = vectors->cords;
        while (tmp != NULL)
        {
            tmp->value = 0;
            tmp = tmp->next;
        }

        vectors = vectors->next;
    }
}

/*add vector 2 into vector 1*/
void add_vectors(struct vector *vec1,struct vector *vec2)
{
    struct cord *cord1 = vec1->cords, *cord2 = vec2->cords;
    while (cord1 != NULL)
    {
        cord1->value += cord2->value;
        cord1 = cord1->next;
        cord2 = cord2->next;
    }
}

/*turn sum of vectors into average*/
void average_vectors(struct vector *centroids_sum,int *centroids_size, struct vector *old_centroids)
{
    int i = 0;
    struct cord *curr_cord, *old_curr_cord;

    while(centroids_sum != NULL)
    {
        curr_cord = centroids_sum->cords;

        if(centroids_size[i] != 0)
            while(curr_cord != NULL)
            {
                curr_cord->value /= (double)centroids_size[i];
                curr_cord = curr_cord->next;
            }
        else{
            old_curr_cord = old_centroids->cords;
            while(curr_cord != NULL)
            {
                curr_cord->value = old_curr_cord->value;
                curr_cord = curr_cord->next;
                old_curr_cord = old_curr_cord->next;
            }
        }
        centroids_sum = centroids_sum->next;
        old_centroids = old_centroids-> next;
        i++;
    }
}


/* free connected cords */
void free_cords(struct cord *cords) {
    struct cord *tmp;
    while (cords != NULL) {
        tmp = cords;
        cords = cords->next;
        free(tmp);
    }
}

/* free connected vectors*/
void free_vector(struct vector *current) {
    struct vector *tmp;
    while (current != NULL) {
        tmp = current;
        current = current->next;
        free_cords(tmp->cords);
        free(tmp);
    }
}


/*initialize centroid and centroidsum from the first K vectors in database - centroid_sum values will be initialized to 0*/
struct vector *initialize_centroids_sum(struct vector *database,int K)
{
    int i;
    struct vector *current = NULL,*centroid_head = NULL, *copy = NULL;

    for(i = 0;i < K;i++)
    {
        if(database == NULL)
        {
            free_vector(centroid_head);
            return NULL;
        }
        copy = copy_vector_from(database);
        if(copy == NULL) /* faild aloocating memory*/
        {
            if(centroid_head != NULL)
                free_vector(centroid_head);
            return NULL;
        }
        copy->next = NULL;
        

        if(centroid_head == NULL)
        {
            centroid_head = copy;
            current = copy;
        } else {
            current->next = copy;
            current = current->next;
        }
        database = database->next;
    }

    return centroid_head;
}

struct vector *copy_vector_from(struct vector *other)
{
    struct vector *tmp;
    struct cord *current_cord, *other_cord;

    tmp = malloc(sizeof(struct vector));
    if(tmp == NULL)
        return NULL;
    tmp->next = NULL;
    tmp->cords = malloc(sizeof(struct cord));
    if(tmp->cords == NULL)
    {
        free(tmp);
        return NULL;
    }
    
    current_cord = tmp->cords;
    other_cord = other->cords;

    while (other_cord != NULL) 
    {
        current_cord->value = other_cord->value;
        if(other_cord->next == NULL)
        {
            current_cord->next = NULL;
        } else {
            current_cord->next = malloc(sizeof(struct cord));
            if(current_cord->next == NULL)
            {
                free_vector(tmp);
                return NULL;
            }
        }
        current_cord = current_cord->next;
        other_cord = other_cord->next;
    } 
    return tmp;
}


/* wrapper function for fit */
static PyObject* fit(PyObject *self, PyObject *args) {
    
    
    PyObject *centroids_table, *data_base_table;
    PyObject *result, *tmp_row, *python_double;
    int K, iter, i, j, n, m;
    double eps;
    struct vector *data_base0, *centroids0, *outcome, *curr_vec; 
    struct cord *curr_cord;

    if (!PyArg_ParseTuple(args, "OOiid",&data_base_table, &centroids_table, &K, &iter, &eps)) {
        return NULL;
    }
    
    
    data_base0 = initialize_data(data_base_table);
    if(data_base0 == NULL)
    {
        return NULL;
    }
        
    centroids0 = initialize_data(centroids_table);
    if(centroids0 == NULL)
    {
        free_vector(data_base0);
        return NULL;
    }
    if(K >= PyObject_Length(data_base_table))
    {
        PyErr_SetString(PyExc_ValueError, "Invalid number of clusters!\n");
        free_vector(centroids0);
        free_vector(data_base0);
        return NULL;
    }
    

    outcome = fit_func(data_base0, centroids0, K, iter, eps);
    if(outcome == NULL)
        return NULL;
        

    n = PyObject_Length(centroids_table);
    m = PyObject_Length(PyList_GetItem(centroids_table, 0));
    result = PyList_New(n);
    curr_vec = outcome;

    /* convert from vectors to PyList*/
    for(i = 0; i < n; i++){
        tmp_row = PyList_New(m);
        curr_cord = curr_vec->cords;
        for(j = 0; j < m; j++){
            python_double = Py_BuildValue("d",curr_cord->value);
            PyList_SetItem(tmp_row, j, python_double);
            curr_cord = curr_cord->next;
        }
        PyList_SetItem(result, i, tmp_row);
        curr_vec = curr_vec->next;
    }
    
    free_vector(outcome);
    return result;
}


struct vector *initialize_data(PyObject *table)
{
    PyObject *lst;
    PyObject *item;

    struct vector *head_vec, *curr_vec;
    struct cord *head_cord, *curr_cord;
    double num;
    int rows_num,cols_num,i,j;

    rows_num = PyObject_Length(table); /* save length of rows*/
    if (rows_num < 0) {
        return NULL;
    }

    lst = PyList_GetItem(table, 0);
    if (lst == NULL) {
        return NULL;
    }
    cols_num = PyObject_Length(lst); /* save length of columns*/
    if (cols_num < 0) {
        return NULL;
    }

    /* create first vector */
    head_vec = malloc(sizeof(struct vector));
    if(head_vec == NULL)
        return NULL;
    curr_vec = head_vec;
    curr_vec->next = NULL;

    for(i = 0; i < rows_num; i++)
    {   /* get next list */
        lst = PyList_GetItem(table, i);
        if (lst == NULL) {
            return NULL;
        }
        /* create current cord chain */
        head_cord = (struct cord *)malloc(sizeof(struct cord));
        if(head_cord == NULL)
        {
            free_vector(head_vec);
            return NULL;
        }
        curr_cord = head_cord;
        curr_cord->next = NULL;

        for(j = 0; j < cols_num; j++)
        {   /* get value */
            item = PyList_GetItem(lst, j);
            num = PyFloat_AsDouble(item);
            /*attach value */
            curr_cord->value = num;
            if (j < cols_num - 1) {
                curr_cord->next = malloc(sizeof(struct cord));
                if(curr_cord->next == NULL)
                {
                    free_vector(head_vec);
                    return NULL;
                }
                curr_cord = curr_cord->next;
                curr_cord->next = NULL;
            } else {
                curr_cord->next = NULL;
            }
        }
        /* attach the created cords */
        curr_vec->cords = head_cord;

        /* alocate new vector if there is more */
        if (i < rows_num -1) {
            curr_vec->next = malloc(sizeof(struct vector));
            if(curr_vec->next == NULL)
            {
                free_vector(head_vec);
                return NULL;
            }
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
        } else { /* last vector */
            curr_vec->next = NULL;
        }
    }

    return head_vec;
}


static PyMethodDef kmeansMethods[] = {
    {"fit",                   /* the Python method name that will be used */
      (PyCFunction) fit, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("fit calculate K centroids of given data\n"
        "it expect to get:"
        "1 - list of lists - each inner list represent a point in the space\n"
        "2 - list of lists - each inner list represent a point of the K centroids\n"
        "3 - integer k - as the k centroids it expects\n"
        "4 - integer iter - as the amount of repeating the code\n"
        "5 - double epsilon - The range of error difference of the centroids\n"
        "its return the k computed cetroids as a list of lists\n")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     
};

static struct PyModuleDef kmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    kmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}

