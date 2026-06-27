import sys
import numpy as np
import pandas as pd
import math # check if to delete
import mykmeanssp as kmeans

# calculate distance in space for 2 vectors -> insert 2 vectors(2 dim array)
def distance (current, other):
    vector = 0.0
    for x,y in zip(current, other):
        vector += (x - y) ** 2
    return vector ** 0.5

def closest(vector, centroids):
    the_closest = 0
    min_length = distance(vector, centroids[0])
    for i in range(1, len(centroids)):
        current_length = distance(vector, centroids[i])
        if min_length > current_length:
            min_length = current_length
            the_closest = i
    return the_closest, min_length

### main ###

# check for legal inputs
if len(sys.argv) != 5 and len(sys.argv) != 6:
    print("An Error Has Occurred")
    sys.exit(1)

try: #K
    if not float(sys.argv[1]).is_integer():
        print("Invalid number of clusters!")
        exit(1)
except ValueError:
    print("Invalid number of clusters!")
    exit(1)
try: #iter
    if len(sys.argv) == 6:
        if not float(sys.argv[2]).is_integer():
            print("Invalid maximum iteration!")
            sys.exit(1)
except ValueError:
    print("Invalid maximum iteration!")
    exit(1)
try: #EPSILON
    if len(sys.argv) == 6:
        float(sys.argv[3])
    else:
        float(sys.argv[2])
except ValueError:
    print("Invalid epsilon!")
    exit(1)

#apply values
if len(sys.argv) == 6:
    K = int(float(sys.argv[1]))
    iter = int(float(sys.argv[2]))
    EPSILON = float(sys.argv[3])
    adress1 = sys.argv[4]
    adress2 = sys.argv[5]
else: # == 5
    K = int(float(sys.argv[1]))
    iter = 300
    EPSILON = float(sys.argv[2])
    adress1 = sys.argv[3]
    adress2 = sys.argv[4]

# read the CSV files
part1 = pd.read_csv(adress1 ,header = None)
part2 = pd.read_csv(adress2 ,header = None)

combined = pd.merge(part1,part2,on = 0).sort_values(by = 0,ascending = True).reset_index(drop = True).drop(0,axis = 1)

# errors in input handeling
if K <= 1 or combined.shape[0] <= K:
    print("Invalid number of clusters!")
    sys.exit(1)
if iter <= 1 or 1000 <= iter:
    print("Invalid maximum iteration!")
    sys.exit(1)
if EPSILON < 0:
    print("Invalid epsilon!")
    sys.exit(1)

np.random.seed(1234)
picked = set()
picked_str = ""

#pick 1 random vector as the first centroid
rnd = np.random.choice(np.arange(combined.shape[0]))
picked_str += str(rnd)
first_centroid = combined.iloc[rnd].copy().tolist()
centroids = pd.DataFrame([first_centroid])

#create the first k centroids
for _ in range(K-1):
    distances = []
    sum_distance = 0.0
    for i in range(combined.shape[0]): # calculate probabilities
        if i in picked:
            distances.append(0)
            continue

        idx, curr_dis = closest(combined.iloc[i].tolist(), centroids.values.tolist())
        distances.append(curr_dis)
        sum_distance += curr_dis
        
    probabilities = [distances[i] / sum_distance for i in range(len(distances))]

    if not math.isclose(sum(probabilities), 1.0): # check for error in calculating the probabilities
        print("summing probabilities isnt 1!!!")
        sys.exit(1)
    
    chosen = np.random.choice(range(combined.shape[0]), size=None, replace=False,p=probabilities) # pick the next centroid
    picked_str += "," + str(chosen)
    picked.add(chosen)
    # take out the new centroid and save him
    deleted = combined.iloc[chosen].copy().tolist()
    centroids = pd.concat([centroids,pd.DataFrame([deleted])],ignore_index=True)
    
try:
    result = kmeans.fit(combined.values.tolist(),centroids.values.tolist(),K,iter,EPSILON)

    print(picked_str)
    # print the centroids
    for vector in result:
        row = ""
        for value in vector:
            row += ("%.4f" % value) + ","
        print(row[:-1])
except ValueError as e:
    print(e)
    sys.exit(1)
except Exception as e:
    print("An Error Has Occurred")
    sys.exit(1)