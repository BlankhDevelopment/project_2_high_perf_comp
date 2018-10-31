#!/bin/bash
#SBATCH --job-name=proj2
#SBATCH --output=res.txt
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:59:00

export OPENBLAS_NUM_THREADS=1
./cs211_proj2_2 1001
./cs211_proj2_2 2001
./cs211_proj2_2 3001
./cs211_proj2_2 4001
./cs211_proj2_2 5001
