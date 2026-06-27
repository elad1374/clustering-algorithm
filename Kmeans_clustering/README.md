# K-Means++ Clustering Optimization

## Overview
This repository contains a high-performance, hybrid **C/Python** implementation of the **K-Means++** clustering algorithm. 

The project strategically splits the computational workload: Python handles data preprocessing, merging datasets, and executing the smart probabilistic initialization (K-Means++). Once the initial centroids are established, the system delegates the heavy iterative distance computations to a custom-built **C extension** for maximum execution speed.

## Project Structure
* `kmeansmodule.c`: The core C extension written from scratch. It implements the standard K-Means iterative assignment and update steps using memory-efficient linked lists (`struct vector` and `struct cord`).
* `kmeanspp.py`: The Python driver that handles argument parsing, data merging (via Pandas), K-Means++ initialization, and interfacing with the C module.
* `setup.py`: The standard setuptools configuration script used to compile the C code into a dynamic Python library (`mykmeanssp`).

## Algorithmic Pipeline
1. **Data Joining:** Reads two distinct CSV datasets, performs an inner merge on their primary keys, and sorts the unified dataset.
2. **K-Means++ Initialization (Python):** * Selects the first centroid uniformly at random.
   * For every subsequent centroid (up to $K$), computes the distance $D(x)$ from each data point to its closest existing centroid.
   * Selects the next centroid probabilistically, where the probability of choosing point $x$ is proportional to its squared distance:
     $$P(x) = \frac{D(x)^2}{\sum_{i} D(x_i)^2}$$
3. **K-Means Optimization (C Engine):** Iteratively assigns points to the nearest centroid using Euclidean distance and recalculates the centroid means until convergence (either maximum iterations are reached or centroid movement is strictly less than $\epsilon$).

## Setup & Compilation

### 1. Prerequisites
Ensure your environment contains the Python C development headers and required data libraries:
```bash
pip install numpy pandas