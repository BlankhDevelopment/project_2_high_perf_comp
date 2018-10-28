#include <stdio.h>                                                                             
#include "cblas.h"                                                                            
#include "lapacke.h"                                                                           
#include <string.h>     
#include <iostream>                                                                       

                                                                                               

using namespace std;                                                                           

                                                                                               

int main (int argc, const char * argv[]) {                                                     

                                                                                               

  int m = 3; //# of rows                                                                                  
  int n = 3; //# of columns                                                                                   
  int lda = 3; //length of first dimenssion                                                                                 
  int ldb = 3;                                                                                 
                                                                                               

  double *A; 
  double *B;
  int *ipiv;

  A = (double*)malloc(m*n * sizeof(double));
  B = (double*)malloc(m * sizeof(double));
  ipiv = (int*)malloc(m * sizeof(int));                                                                     

                                                                                                                                                                                                                                                                                                                                                                                           

  A[0] = 4; A[1] = 6;  A[2] =  7;        B[0] = 2; B[1] = 10; B[2] = 6;                                                          
  A[3] = 3; A[4] = 3;  A[5] = 4;         B[3] = 3; B[4] = 5; B[5] = 4;                               
  A[6] = 6; A[7] = 4;  A[8] = 6;         B[6] = 11; B[7] = 6; B[8] = 11;

  //dtrsm paramters
  CBLAS_SIDE SIDE = 'l';
  CBLAS_UPLO UPLO = 'u';
  CBLAS_TRANSPOSE TRANSA = 'n';
  CBLAS_DIAG DIAG = 'u';
  double ALPHA = 0;                                                                                                                                                                                                                                                                                    

  LAPACKE_dgetrf( LAPACK_COL_MAJOR, m, n, A, lda, ipiv );  
  cblas_dtrsm( LAPACK_COL_MAJOR, SIDE, UPLO, TRANSA, DIAG, m, n, ALPHA, A, lda, B, ldb);                                    

                                                                                               

  for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }                                                                                            
   
   cout << endl << endl;
  
    cblas_dtrsm(LAPACK_COL_MAJOR, SIDE, UPLO, TRANSA, DIAG, m, n, ALPHA, A, lda, B, ldb);  

   for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }       

    cout << endl << endl;
  
     cblas_dtrsm(LAPACK_COL_MAJOR, SIDE, UPLO, TRANSA, DIAG, m, n, ALPHA, A, lda, B, ldb);      

   for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }       
  
                                                                                               
  /*cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);      
for(int i=0; i < n; i++){
     printf("  %lf ", B[i]);
     printf("\n");
}

*/
  return 0;

}
