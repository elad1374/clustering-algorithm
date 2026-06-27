import sys
import numpy as np
import pandas as pd
import math 
import symnmf


# constant
ERROR = "An Error Has Occurred"

# printing 2 dim lists
def print_matrix(vectors):
    for vector in vectors:
        row = ""
        for value in vector:
            row += "%.4f" % value + ","
        print(row[:-1])

# program failed - prints ERROR massage and exit the program
def failed():
    print(ERROR)
    exit(1)


### Main ###
try: 
    #read inputs
    if len(sys.argv) == 4: # K included
        K = int(float(sys.argv[1]))
        goal = sys.argv[2]
        vectors_list = pd.read_csv(sys.argv[3], header = None).values.tolist()
    else: # K not included
        goal = sys.argv[1]
        vectors_list = pd.read_csv(sys.argv[2], header = None).values.tolist()

    if(goal != "sym" and goal !="ddg" and goal != "norm" and goal != "symnmf"): # ensure goal viability
        failed()

    # perform goal
    if goal == "sym": # create the similarity matrix A
        matrix = symnmf.sym(vectors_list)
        print_matrix(matrix)

    if goal == "ddg": # create the diagonal degree matrix D
        matrix = symnmf.ddg(vectors_list)
        print_matrix(matrix)

    if goal == "norm": # create the normalized matrix W
        matrix = symnmf.norm(vectors_list)
        print_matrix(matrix)

    if goal == "symnmf": # perform all symnmf algorithem
        # create the normalized W matrix
        matrix_W = symnmf.norm(vectors_list) 
        # create the randomized H matrix
        np.random.seed(1234)
        high_boundry = 2 * math.sqrt((np.array(matrix_W).mean()) / K)
        matrix_H = np.random.uniform(low = 0.0, high = high_boundry, size = (len(vectors_list), K)).tolist()
        # create the clustering matrix H and print
        new_matrix_H = symnmf.symnmf(matrix_H, matrix_W) 
        print_matrix(new_matrix_H)
        
except: # an error occurred, print and terminate
    failed()
