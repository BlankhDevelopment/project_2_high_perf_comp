#!/bin/bash
#SBATCH --job-name=proj2
#SBATCH --output=res.txt
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:59:00

export OPENBLAS_NUM_THREADS=1
./cs211_proj2_1_1 1000
./cs211_proj2_1_1 2000
./cs211_proj2_1_1 3000
./cs211_proj2_1_1 4000
./cs211_proj2_1_1 5000
