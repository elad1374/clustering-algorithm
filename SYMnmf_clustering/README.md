# Symmetric Non-Negative Matrix Factorization (SymNMF) Clustering

## Overview
This repository contains a hybrid **C/Python** implementation of the **Symmetric Non-negative Matrix Factorization (SymNMF)** algorithm, alongside a comprehensive performance benchmark against the classic **K-Means** clustering algorithm using **Silhouette Scores**.

SymNMF is an advanced unsupervised learning technique used for graph clustering and community detection. It factorizes a symmetric similarity matrix derived from data points into a low-rank, non-negative cluster-assignment matrix *H*, yielding soft clustering coordinates that are subsequently discretized.

## Project Structure
The project combines the raw computational efficiency of **C** for matrix calculations with the ease of analysis provided by **Python**:

* `symnmf.c` / `prototypes.h`: The core C engine containing execution routines for similarity modeling, diagonal degree computations, matrix normalization, and the iterative multiplicative update rules for $H$.
* `symnmfmodule.c`: The CPython extension bridge wrapping native C functions into a Python-callable binary module (`symnmf`).
* `symnmf.py`: The command-line driver executing specific operational targets (`sym`, `ddg`, `norm`, `symnmf`) using Pandas and NumPy.
* `analysis.py`: The benchmarking framework that fits both SymNMF and K-Means onto an input dataset, evaluates cluster qualities, and prints their respective **Silhouette Scores** via `scikit-learn`.
* `setup.py`: The standard setuptools configuration compilation script for compiling the C extension wrapper.
* `Makefile`: Builds the standalone C executable natively.

## Mathematical Pipeline
Given a set of $N$ points in a $d$-dimensional space, the algorithm executes the following structural phases:

1. **Similarity Matrix ($A$):** Computes pairwise local Euclidean distances using a Gaussian kernel:
   $$A_{ij} = \exp\left(-\frac{||x_i - x_j||^2}{2}\right) \quad \text{for } i \neq j, \quad A_{ii} = 0$$
2. **Diagonal Degree Matrix ($D$):** Constructs a diagonal matrix where each $D_{ii}$ represents the sum of row $i$ in matrix $A$:
   $$D_{ii} = \sum_{j} A_{ij}$$
3. **Normalized Similarity Matrix ($W$):** Computes a symmetrically normalized graph Laplacian layout:
   $$W = D^{-1/2} A D^{-1/2}$$
4. **Factorization Optimization ($H$):** Iteratively solves the objective function by updating $H$ using Multiplicative Update Rules until convergence ($\le 300$ iterations or within $\epsilon = 10^{-4}$):
   $$\min_{H \ge 0} ||W - HH^T||_F^2$$

## Setup & Compilation

### 1. Prerequisites
Ensure your environment contains the Python C development headers, standard build utilities, and key data libraries:
```bash
pip install numpy pandas scikit-learn

# Build the Python C extension
python setup.py build_ext --inplace

# Build the C executable
make

# Output the Similarity Matrix (A)
python symnmf.py sym sample_data.txt

# Output the Diagonal Degree Matrix (D)
python symnmf.py ddg sample_data.txt

# Output the Normalized Similarity Matrix (W)
python symnmf.py norm sample_data.txt

# Full SymNMF Clustering
python symnmf.py 3 symnmf sample_data.txt

# Performance Analysis
python analysis.py 3 sample_data.txt

