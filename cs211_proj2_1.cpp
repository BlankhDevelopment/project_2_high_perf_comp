#include <stdio.h>                                                                             

#include "cblas.h"                                                                             

#include "lapacke.h"                                                                           

#include <string.h>                                                                            

                                                                                               

using namespace std;                                                                           

                                                                                               

int main (int argc, const char * argv[]) {                                                     

                                                                                               

  int m = 2;                                                                                   

  int n = 2;                                                                                   

  int lda = 2;                                                                                 

  int ldb = 2;                                                                                 

                                                                                               

  double * A = new double[m*n];                                                                

  double * B = new double[m];                                                                  

  int * ipiv = new int[m];                                                                     

                                                                                                                                                                                                                                                                                                                                                                                           

  A[0] = 4; A[1] = 6;  B[0] = 2;                                                               

  A[2] = 3; A[3] = 3;  B[1] = 10;                                                              

                                                                                                                                                                                                                                                                                         

  LAPACKE_dgetrf( LAPACK_COL_MAJOR, m, n, A, lda, ipiv );                                      

                                                                                               

  for (int i = 0; i < m; i++) {                                                                

      for (int j = 0; j < n; j++){                                                             

          printf("  %lf ", A[lda*j+i]);                                                        

      }                                                                                        

      printf("\n");                                                                            

                                                                                               

  }                                                                                            

                                                                                               

  cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);      

for(int i=0; i < n; i++){

     printf("  %lf ", B[i]);

     printf("\n");

}

  return 0;

}
