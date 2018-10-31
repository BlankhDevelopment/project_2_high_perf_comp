#include <stdio.h>         
#include <ctime>
#include <iomanip>
#include <cstdlib>               
#include <math.h>                                                     
#include "cblas.h"                                                                            
#include "lapacke.h"                                                                           
#include <string.h>     
#include <iostream>                                                                                                                                                       

using namespace std; 


//THIS IS A CAREFUL CURATION OF THE MATLAB CODE

int* mydegtrf(double arr[], int array_size)
{
    int block = 10;
    int n = array_size;
    int *pvt;
    pvt = (int*)malloc((array_size*array_size) * sizeof(int));  
    for(int i = 0; i < array_size; i++)
    {
        pvt[i] = i;
    }
    double *tempv;
    tempv = (double*)malloc(array_size* sizeof(double)); 
    for(int i = 0; i < n-1; i+=block) //for 1-->(n-1)
    {
        int maxind = i;
        double max = fabs(arr[i*n + i]); //this will be the same as A(i,i)
        for (int t = i + 1; t < n; t += block)
        {
            for (int i1 = i; i1 < i + block; i1++)
            {
                for (int t1 = t; t1 < t + block; t1++)
                {
                    if(fabs(arr[t1*n + i1])>max)
                    {
                        maxind = t1; 
                        max = fabs(arr[t1*n+i1]);
                    }
                }
            }
        }
        if (max == 0)
        {
             cout << "LUFactoration failed: coefficient matrix is singular" << endl; return 0;
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
                for(int r = 0; r < array_size; r+=block)
                {
                    for(int i1 = i; i1 < i + block; i1++)
                    {
                        for(int r1 = r; r1 < r + block; r1++)
                        {
                            tempv[r1] = arr[i1*n+r1];
                            arr[i1*n + r1] = arr[maxind*n + r1];
                            arr[maxind*n + r1] = tempv[r1];
                        }
                    }
                }
            }
        }
        
       //blocked code
        for (int j = i + 1; j < n; j += block)
        {
            for (int k = i + 1; k < n; k += block)
            {
                for(int i1 = i; i1 < i + block; i1++)
                {
                    for(int j1 = j; j1 < j + block; j1++)
                    {
                        register double r2 = j1*n + i1;
                        register double r3 = i1*n + i1;
                        arr[r2] = arr[r2]/arr[r3];
                        register double r1 = arr[j1*n + i1];
                        for(int k1 = k; k1 < k + block; k1++)
                        {
                            register int r4 = j1*n + k1;
                            register int r5 = i1*n + k1;
                            register double r6 = r1 * arr[r5];
                            register double r7 = arr[r4];
                            arr[r4] = r7 - r6;
                        }
                    }
                }
            }
        }
        
        /*
        //factorization
        for(int j = i + 1; j < n; j++)
        {
            arr[j*n + i] = arr[j*n + i]/arr[i*n + i];
            for(int k = i + 1; k < n; k++)
            {
                arr[j*n + k] = arr[j*n + k] - (arr[j*n+i] * arr[i*n + k]);
            }
        }
        */
        
    }

    return pvt;
}

//FORWARD SUBSTITUTION curated with Matlab code
double* mydtrsmfwd(double arr[], double arr2[], int pivot[], int array_size)
{
    double sum = 0.0;
    int n = array_size;
    double *y;
    y = (double*)malloc(array_size* sizeof(double));
    y[0] = arr2[pivot[0]];
    for (int i = 1; i < n; i++)
    {
        sum = 0.0;
        for(int r = 0; r < i; r++)
        {
            sum = sum + (y[r] * arr[i*n + r]);
        }
        y[i] = arr2[pivot[i]] - sum;
    }
    return y;
}

//BACKWARDS SUBSTITUTION curated with Matlab code
double* mydtrsmbwd(double arr[], double arr2[], int array_size)
{
    double sum = 0.0;
    int n = array_size;
    double *x;
    x = (double*)malloc(array_size* sizeof(double));
    x[n-1] = arr2[n-1] / arr[(n-1)*(n) + (n-1)]; //x(n) = y(n) / A(n,n)
    /*

    */
    for (int i = n-2; i >= 0; i--)
    {
        sum = 0.0;
        for(int r = i + 1; r < n; r++)
        {
            sum += x[r] * arr[i*n+r];
        }
        x[i] = (arr2[i] - sum) / arr[i*n + i];
    }
    return x;
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
  int ldb = 1; //leading dimension of b                                                                    
                                                                                               
  double *A; 
  double *B;
  double *C;
  int *ipiv;
  
  A = (double*)malloc(m*n * sizeof(double));
  B = (double*)malloc(m * sizeof(double));
  C = (double*)malloc(m * sizeof(double));
  ipiv = (int*)malloc(m* sizeof(int));             

  FillMatrix(A, B, n);     
/*
  A[0] = 1; A[1] = 1; A[2] = 1; A[3] = 4;
  A[4] = 3; A[5] = -1; A[6] = 3; A[7] = 5;
  A[8] = 3; 
                                      
 B[0] = 1; B[1] = 6; B[2] = 4;
// C[0] = 6; C[1] = -0.5; C[2] = -0.454545;
*/
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

cout << "Performing blocked mydgetrf (LU factorization) with matrix size: " << n << endl;

clock_t t;

t = clock();
ipiv = mydegtrf(A, n);
t = clock() - t;

cout << "This process took: " << (double(t) / CLOCKS_PER_SEC) << " seconds" << endl;


//outputs matrix A
/*
 for (int i = 0; i < m; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", A[lda*j+i]);                                                        
      }                                                                                        
      printf("\n");                                                                                                                                                                   
  }    
*/
/*  //swapping arrays 
for(int i = 0; i < n; i++)
{
    double temp = B[i];
    B[i] = B[ipiv[i]];
    B[ipiv[i]] = temp;
}
*/ 
/*
cout << "PIVOT ARRAY: " << endl;
for (int i = 0; i < n; i++)
{
    cout << ipiv[i] << endl;
}
//ipiv[2] = 2;
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
    
    
   
     cout << "Performing mydtrsm (fwd/bwd substitution) with matrix size: " << n << endl;
    
      t = clock();
     C = mydtrsmfwd(A, B, ipiv, n); //performing forward subsitution
  /*  cout << endl << "AFTER PERFORMING FORWARD SUBSTITUTION ON MATRIX B, WE GET: " << endl;                                                           
      for (int j = 0; j < n; j++)
      {                                                             
          cout << C[j];
          cout << endl;                                                       
      }     
      */  
      C = mydtrsmbwd(A, C, n);  // performing backward substitution

      t = clock() - t;

      cout << "This process took: "  << (double(t) / CLOCKS_PER_SEC) << " seconds" << endl;

/*
        cout << endl << "AFTER PERFORMING BACKWARDS SUBSTITUTION ON MATRIX B, WE GET: " << endl;                                                           
      for (int j = 0; j < n; j++)
      {                                                             
          cout << C[j];
          cout << endl;                                                       
      }                                                                                        
      */                                                                                                                                                               
   

/*

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
