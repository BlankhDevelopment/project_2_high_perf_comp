#include <stdio.h>                                                                             
#include "cblas.h"                                                                            
#include "lapacke.h"                                                                           
#include <string.h>     
#include <iostream>                                                                       

                                                                                               

using namespace std;                                                                           

                                                                                               

int main (int argc, const char * argv[]) {                                                     

                                                                                               

  int m = 4; //# of rows                                                                                  
  int n = 4; //# of columns                                                                                   
  int lda = 4; //length of first dimenssion                                                                                 
  int ldb = 4;                                                                                 
                                                                                               

  double *A; 
  double *B;
  double *C;
  int *ipiv;

  A = (double*)malloc(m*n * sizeof(double));
  B = (double*)malloc(m*n * sizeof(double));
  C = (double*)malloc(m*n * sizeof(double));
  ipiv = (int*)malloc(m* sizeof(int));                                                                     

  A[0] = 2; A[1] = 0; A[2] = 0; A[3] = 0;
  A[4] = -2; A[5] = 1; A[6] = -1; A[7] = 0;
  A[8] = 6; A[9] = 2; A[10] = 1; A[11] = 0;
  A[12] = 0; A[13] = 0; A[14] = 0; A[15] = 0;                                                                                                                                                                                                                                                                                                                                                                                       

 B[0] = 4; B[1] = -4; B[2] = 15; B[3] = -1;                                                         


cout << endl << "OUTPUTTING MATRIX A: " << endl; 


 for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }                 
  cout << endl << endl;
    
      LAPACKE_dgetrf( LAPACK_COL_MAJOR, m, n, A, lda, ipiv );  
     // cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, 1.0, A, lda, B, ldb);   

cout << endl << "AFTER PERFORMING LU FACTORIZATION ON MATRIX A, WE GET: " << endl;                                                                 

  for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }                                                                                            
   
   cout << endl << endl;

  
      cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, 1.0, A, lda, B, ldb);    

                                                                 
      for (int j = 0; j < n; j++)
      {                                                             
          printf("  %lf ", B[lda*j+i]); 
          cout << endl;                                                       
      }                                                                                        
                                                                                                                                                                      
    

    cout << endl << endl;

  
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, m, n, 1.0, A, lda, B, ldb);     

                                                              
      for (int j = 0; j < n; j++)
      {                                                             
          printf("  %lf ", B[lda*j+i]);                                                        
          cout << endl;
      }                                                                                   
                                                                                                                                                 
        
  
                                                                                               
  /*cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);      
for(int i=0; i < n; i++){
     printf("  %lf ", B[i]);
     printf("\n");
}

*/
  return 0;
}
