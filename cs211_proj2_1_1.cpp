#include <stdio.h>         
#include <ctime>
#include <iomanip>
#include <cstdlib>                                                                    
#include "cblas.h"                                                                            
#include "lapacke.h"                                                                           
#include <string.h>     
#include <iostream>                                                                       

                                                                                               

using namespace std;                                                                           

void FillMatrix(double arr1[], double arr2[], int array_size)
{
	int m = 0;

     srand(1);
	for (m = 0; m < array_size*array_size; ++m)
	{
		arr1[m] = rand() % 100;
	}

    for (m = 0; m < array_size; m++)
    {
        arr2[m] = rand() % 100;
    }
}
                                                                                               

int main (int argc, const char * argv[]) {                                                     
                                                                                            
  const int n = atoi(argv[1]);
  int m = n; //# of rows                                                                                                                                                                  
  int lda = n; //length of first dimenssion                                                                                 
  int ldb = n;                                                                                 
                                                                                               
  double *A; 
  double *B;
  double *C;
  int *ipiv;

  A = (double*)malloc(m*n * sizeof(double));
  B = (double*)malloc(m * sizeof(double));
  C = (double*)malloc(m * sizeof(double));
  ipiv = (int*)malloc(m* sizeof(int));             

  FillMatrix(A, B, n);     


  A[0] = 4; A[1] = 6; A[2] = 3; A[3] = 3;
 /* A[4] = -2; A[5] = 1; A[6] = -1; A[7] = 0;
  A[8] = 6; A[9] = 2; A[10] = 1; A[11] = 0;
  A[12] = 0; A[13] = 0; A[14] = 0; A[15] = 0;                                                                                                                                                                                                                                                                                                                                                                                       
 B[0] = 4; B[1] = -4; B[2] = 15; B[3] = -1;    */                                         
 B[0] = 2; B[1] = 10;
 C[0] = 2; C[1] = 10;


/*
cout << endl << "OUTPUTTING MATRIX A: " << endl; 


 for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }                 
  cout << endl;
*/

cout << "Performing dgetrf (LU factorization) with matrix size: " << n << endl;

clock_t t;

t = clock();
    
LAPACKE_dgetrf( LAPACK_COL_MAJOR, m, n, A, lda, ipiv );  

t = clock() - t;

cout << "This process took: " << (double(t) / CLOCKS_PER_SEC) << " seconds" << endl;




  cout << endl << "AFTER PERFORMING LU FACTORIZATION ON MATRIX A, WE GET: " << endl;                                                                 

  for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }                                                                                            
     
    
    
    
     cout << "Performing dtrsm (fwd/bwd substitution) with matrix size: " << n << endl;
    
      t = clock();
  
      cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);    


        cout << endl << "AFTER PERFORMING FORWARD SUBSTITUTION ON MATRIX B, WE GET Y-coefficient list: " << endl;                                                           
      for (int j = 0; j < n; j++)
      {                                                             
          cout << B[j];
          cout << endl;                                                       
      }                                                                                        
      
      B[0] = C[0];
      B[1] = C[1];


    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);     

    
      t = clock() - t;

      cout << "This process took: "  << (double(t) / CLOCKS_PER_SEC) << " seconds" << endl;

     cout << endl << "AFTER PERFORMING BACKWARD SUBSTITUTION ON MATRIX B, WE GET X-coefficient list: " << endl;
                                                            
      for (int j = 0; j < n; j++)
      {                                                             
          cout << B[j];                                                      
          cout << endl;
      }                                                                                   
      

  return 0;
}
