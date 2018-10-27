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
  int ldb = 2;                                                                                 
                                                                                               

  double *A; 
  double *B;
  int *ipiv;

  A = (double*)malloc(m*n * sizeof(double));
  B = (double*)malloc(m * sizeof(double));
  ipiv = (int*)malloc(m * sizeof(int));                                                                     

                                                                                                                                                                                                                                                                                                                                                                                           

  A[0] = 4; A[1] = 6;  A[2] =  7;        B[0] = 2;                                                               
  A[3] = 3; A[4] = 3;  A[5] = 4;         B[1] = 10;                                                              
  A[6] = 6; A[7] = 4;  A[8] = 6;

  //dtrsm paramters
  char L = 'l';
  char U = 'u';
  char N = 'n';
                                                                                                                                                                                                                                                                                       

  LAPACKE_dgetrf( LAPACK_COL_MAJOR, m, n, A, lda, ipiv );  
  LAPACKE_dtrsm( L, U, N, U, m, n, A, lda, ipiv );                                    

                                                                                               

  for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }                                                                                            
   
   cout << endl << endl;

    LAPACKE_dtrsm( L, U, N, U, m, n, A, lda, ipiv );    

   for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }       

    cout << endl << endl;
  
    LAPACKE_dtrsm( L, U, N, U, m, n, A, lda, ipiv );    

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
