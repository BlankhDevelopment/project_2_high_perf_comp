#include <stdio.h>                                                                             
#include "cblas.h"                                                                            
#include "lapacke.h"                                                                           
#include <string.h>     
#include <iostream>                                                                       

                                                                                               

using namespace std;                                                                           

                                                                                               

int main (int argc, const char * argv[]) {                                                     

                                                                                               

  int m = 2; //# of rows                                                                                  
  int n = 2; //# of columns                                                                                   
  int lda = 2; //length of first dimenssion                                                                                 
  int ldb = 2;                                                                                 
                                                                                               

  double *A; 
  double *B;
  double *C;
  int *ipiv;

  A = (double*)malloc(m*n * sizeof(double));
  B = (double*)malloc(m*n * sizeof(double));
  C = (double*)malloc(m*n * sizeof(double));
  ipiv = (int*)malloc(m* sizeof(int));                                                                     

  A[0] = 4; A[1] = 6; A[2] = 3; A[3] = 3;
/*  A[4] = -2; A[5] = 1; A[6] = -1; A[7] = 0;
  A[8] = 6; A[9] = 2; A[10] = 1; A[11] = 0;
  A[12] = 0; A[13] = 0; A[14] = 0; A[15] = 0;                                                                                                                                                                                                                                                                                                                                                                                       

 B[0] = 4; B[1] = -4; B[2] = 15; B[3] = -1;                                                         
*/
 B[0] = 2; B[1] = 10;

cout << endl << "OUTPUTTING MATRIX A: " << endl; 


 for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }                 
  cout << endl;
    
      LAPACKE_dgetrf( LAPACK_COL_MAJOR, m, n, A, lda, ipiv );  
     // cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, 1.0, A, lda, B, ldb);   

cout << endl << "OUTPUTTING MATRIX B: " << endl; 

 for (int j = 0; j < n; j++)
      {                                                             
          cout << B[j];
          cout << endl;                                                       
      }  

cout << endl;


cout << endl << "AFTER PERFORMING LU FACTORIZATION ON MATRIX A, WE GET: " << endl;                                                                 

  for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }                                                                                            
   
     cout << endl << "AFTER PERFORMING FORWARD SUBSTITUTION ON MATRIX B, WE GET: " << endl;

  
      cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);    

                                                                 
      for (int j = 0; j < n; j++)
      {                                                             
          cout << B[j];
          cout << endl;                                                       
      }                                                                                        
                                                                                                                                                                      
    


    cout << endl << "AFTER PERFORMING BACKWARD SUBSTITUTION ON MATRIX B, WE GET: " << endl;
  
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);     

                                                              
      for (int j = 0; j < n; j++)
      {                                                             
          cout << B[j];                                                      
          cout << endl;
      }                                                                                   
                                                                                                                                                 
  return 0;
}
