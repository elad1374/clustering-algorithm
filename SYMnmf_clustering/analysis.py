import sys
import numpy as np
import pandas as pd
import math 
import symnmf
from sklearn.metrics import silhouette_score

# constants
ERROR = "An Error Has Occurred"
EPSILON = 1e-4
MAX_ITER = 300

# calculate distance in space between current and other vector -> insert 2 vectors(arrays)
def distance(current, other):
    
    vector = 0.0
    for x,y in zip(current, other):
        vector += (x - y) ** 2

    return vector ** 0.5

# return the closest centroid's group index of kmeans
# vector is array, centroids is 2 dim array as the centroids
def closest_centroid_kmeans(vector, centroids):
    
    the_closest = 0
    min_length = distance(vector,centroids[0])
    for i in range(1, K):
        current_length = distance(vector,centroids[i])
        if min_length > current_length:
            min_length = current_length
            the_closest = i
            
    return the_closest

# perform the kmeans algorithem (version hw1) and return the centroids
# vector_list is 2 dim array as the list of vectors(arrays)
def kmeans(vectors_list, K):
    iter = MAX_ITER
    vector_length = len(vectors_list[0])
    vectors_amount = len(vectors_list)

    #initialize the centroids
    centroids = [vectors_list[i].copy() for i in range(K)] # 3 dimenstion array
    clusters_sums = [[0.0 for _ in range(vector_length)] for _ in range(K)] 
    clusters_size = [0 for _ in range(K)] 

    for _ in range(iter):
        
        for i in range(K): # create current empty clusters
            for j in range(vector_length):
                clusters_sums[i][j] = 0
            clusters_size[i] = 0

        
        for vector in vectors_list: # add all vectors to their designated centroids
            idx = closest_centroid_kmeans(vector, centroids) # get the index of the closest centroid
            for i in range(vector_length):
                clusters_sums[idx][i] += vector[i]
            clusters_size[idx] +=1

        
        centroids_too_close = True
        for i in range(K): # compute all new centroids
            if(clusters_size[i] != 0): # not an empty cluster 
                old = centroids[i].copy()
                for j in range(vector_length):
                    centroids[i][j] = clusters_sums[i][j] / clusters_size[i]
                if(distance(old, centroids[i]) > EPSILON):
                    centroids_too_close = False
                
        # check if new centroids are far enough 
        if centroids_too_close:
            break

    return centroids

# get max value's index in array
# returns on empty list -1
def max_index(list):
    try:
        max_idx = 0
        while True:
            if list[max_idx] == max(list):
                break
            max_idx += 1

        return max_idx
    
    except:
        return -1


### Main ###
try: 
    #read inputs    
    K = int(float(sys.argv[1]))
    vectors_list = pd.read_csv(sys.argv[2], header = None).values.tolist()

    # perform symnmf and clustering
    matrix_W = symnmf.norm(vectors_list) 
    np.random.seed(1234)
    high_boundry = 2 * math.sqrt((np.array(matrix_W).mean()) / K)
    matrix_H = np.random.uniform(low = 0.0, high = high_boundry, size = (len(vectors_list), K)).tolist()
    new_matrix_H = symnmf.symnmf(matrix_H, matrix_W) # create the clustering matrix H
    symnmf_clusters = [max_index(new_matrix_H[i]) for i in range(len(vectors_list))]

    # perform kmeans and clustering
    centroids = kmeans(vectors_list, K) # get kmeans centroids
    kmeans_clusters = [closest_centroid_kmeans(vectors_list[i],centroids) for i in range(len(vectors_list))]

    # calculate silhouette scores
    silhouette_score_symnmf = silhouette_score(vectors_list, symnmf_clusters)
    silhouette_score_kmeans = silhouette_score(vectors_list, kmeans_clusters)

    #print results
    print("nmf: %.4f" % silhouette_score_symnmf)
    print("kmeans: %.4f" % silhouette_score_kmeans)
    

except: # an error occurred, print and terminate
    print(ERROR)
    exit(1)