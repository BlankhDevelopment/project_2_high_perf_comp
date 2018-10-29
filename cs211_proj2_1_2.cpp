#include <stdio.h>         
#include <ctime>
#include <iomanip>
#include <cstdlib>                                                                    
#include "cblas.h"                                                                            
#include "lapacke.h"                                                                           
#include <string.h>     
#include <iostream>                                                                       

                                                                                               

using namespace std;           

void mydegtrf(double arr[], int array_size)
{
    int n = array_size;
    int *pvt;
    pvt = (int*)malloc((array_size*array_size) * sizeof(int));  
    double *tempv;
    tempv = (double*)malloc(array_size* sizeof(double)); 
    for(int i = 0; i < n; i++) //for 1-->(n-1)
    {
        int maxind = i;
        int max = abs(arr[i*n + i]); //this will be the same as A(i,i)
        for (int t = i + 1; t <= n; t++)
        {
            if(abs(arr[t*n + i])>max)
            {
                maxind = t; 
                max = abs(arr[t*n+i]);
            }
        
        if (max == 0)
        {
             cout << "LUFactoration failed: coefficient matrix is singular" << endl; return;
        }
        else
        {
            if(maxind != i)
            {
                int temps = pvt[i];
                pvt[i] = pvt[maxind];
                pvt[maxind] = temps;
                //swap rows
                // tempv = A(i,:);
                for(int r = 0; r < array_size; ++r)
                {
                    tempv[r] = A(i*n+r);
                    arr(i*n + r) = arr(maxind*n+r);
                }
            }
        }
        //factorization
        for(int j = i + 1; j <= n; j++)
        {
            arr[j*n + i] = arr[j*n + i]/arr[i*n + i];
            for(int k = i + 1; i <= n; i++)
            {
                arr[j*n + k] = arr[j*n + k] - (arr[j*n+i] * arr[i*n + k]);
            }
        }
    }

}

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
  C = (double*)malloc(m*n * sizeof(double));
  ipiv = (int*)malloc(m* sizeof(int));             

 // FillMatrix(A, B, n);     


  A[0] = 4; A[1] = 6; A[2] = 3; A[3] = 3;
/*  A[4] = -2; A[5] = 1; A[6] = -1; A[7] = 0;
  A[8] = 6; A[9] = 2; A[10] = 1; A[11] = 0;
  A[12] = 0; A[13] = 0; A[14] = 0; A[15] = 0;                                                                                                                                                                                                                                                                                                                                                                                       
 B[0] = 4; B[1] = -4; B[2] = 15; B[3] = -1;  */                                                  
 B[0] = 2; B[1] = 10;


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

cout << "Performing mydgetrf (LU factorization) with matrix size: " << n << endl;

clock_t t;

t = clock();
mydegtrf(A, n);
t = clock() - t;

cout << "This process took: " << (double(t) / CLOCKS_PER_SEC) << " seconds" << endl;

/* outputs matrix B
 for (int j = 0; j < n; j++)
      {                                                             
          cout << B[j];
          cout << endl;                                                       
      }  

cout << endl;
*/

/*
  cout << endl << "AFTER PERFORMING LU FACTORIZATION ON MATRIX A, WE GET: " << endl;                                                                 

  for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }                                                                                            
   */
     
    
    
/*    
     cout << "Performing dtrsm (fwd/bwd substitution) with matrix size: " << n << endl;
    
      t = clock();
  
      cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);    

/*
        cout << endl << "AFTER PERFORMING FORWARD SUBSTITUTION ON MATRIX B, WE GET: " << endl;                                                           
      for (int j = 0; j < n; j++)
      {                                                             
          cout << B[j];
          cout << endl;                                                       
      }                                                                                        
                                                                                                                                                                      
 */   

/*
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);     

    
      t = clock() - t;

      cout << "This process took: "  << (double(t) / CLOCKS_PER_SEC) << " seconds" << endl;

    /* cout << endl << "AFTER PERFORMING BACKWARD SUBSTITUTION ON MATRIX B, WE GET: " << endl;
                                                            
      for (int j = 0; j < n; j++)
      {                                                             
          cout << B[j];                                                      
          cout << endl;
      }                                                                                   
      */

  return 0;
}
