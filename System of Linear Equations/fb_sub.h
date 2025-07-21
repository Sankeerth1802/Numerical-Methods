#include <stdio.h>
#include <stdlib.h>
#define FB_SUB_H

void *forward_sub(float *L,float *C,int n)                             // function that returns a pointer to an array and the arguments are 2 pointers to arrays and an integer(number of columns)
	{
    float *Sol = (float *)malloc(n*sizeof(float));                     // dynamic memory allocation
		
		for (int i=0; i<n; i++)
			{
				*(Sol+i) = 0;                                                  // initialising all elements of solution matrix to zero.
			}
			
		for (int i=0; i<n; i++)                                            // finding x_i = (b_i - summation(j = 1 to i - 1){L[i][j] * x_j}) / L[i][i] where i iterates from 1 to n.
			{
				float sum = 0;
				
				for (int j=0; j<i; j++)
					{
						sum += *(L+(n*i)+j) * *(Sol+j);
					}
				
				*(Sol+i) = (*(C+i) - sum) / *(L+(n*i)+i);
			}
		
		return Sol;                                                        // returning the pointer of solution matrix.
	}

void *backward_sub(float *U,float *C,int n)
	{
		float *Sol = (float *)malloc(n*sizeof(float));
		
		for (int i=0; i<n; i++)
			{
				*(Sol+i)=0;
			}
			
		for (int i=n-1; i>=0; i--)                                         // finding x_i = (b_i - summation(j = i + 1 to n){U[i][j] * x_j}) / U[i][i] but here i iterates from n to 1.
			{
				float sum=0;
				
				for (int j=i+1; j<n; j++)
					{
						sum += *(U+(n*i)+j) * *(Sol+j);
					}
				
				*(Sol+i) = (*(C+i) - sum) / *(U+(n*i)+i);
			}
			
			return Sol;
	}



