cs211_proj2: cs211_proj2_1.cpp
  g++ -o cs211_proj2_1 cs211_proj2_1.cpp -I./lapack-3.8.0/CBLAS/include -I./lapack-3.8.0/LAPACKE/include  ./lapack-3.8.0/liblapacke.a ./lapack-3.8.0/liblapack.a ./lapack-3.8.0/libcblas.a /act/lib/openblas-0.2.20/lib/libopenblas.a -lgfortran -pthread
