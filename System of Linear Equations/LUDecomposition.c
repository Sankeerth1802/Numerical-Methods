#include <stdio.h>
#include "fb_sub.h"

int main()
	{
		FILE *Input,*Output;                                     // file pointers of input and output files.
    
		int n;                                                   // variable to store number of variables.
		
		Input = fopen("InputFile.txt","r");                      // opening input file. 
		
		if (Input == NULL)                                       // error management.
			{
				puts("ERROR: In opening InputFile.txt");
				return 0;
			}
			
				fscanf(Input,"%d",&n);                               // reading the value of number of variables from input file.
				
				float A[n][n],B[n][1];                               // declaring the matrices.
				
				for (int i=0; i<n; i++)
					{
						for (int j=0; j<n; j++)
							{
								fscanf(Input,"%f",&A[i][j]);                 // reading the values from input file.
							}
					}
				
				for (int i=0; i<n; i++)
					{
						fscanf(Input,"%f",&B[i][0]);
					}
				
				fclose(Input);
				
			float L[n][n],U[n][n];                                 // declaring lower and upper triangular matrices.
			
			for (int i=0; i<n; i++)
				{
					L[i][i] = 1;                                       // setting the diagonals of lower triangular matrix as 1.
				}
			
			for (int j=0; j<n; j++)
				{
					for (int i=0; i<j+1; i++)
						{
							if (i != j)
								{
									L[i][j] = 0;                               // setting the elements above the diagonal in lower triangular matrix as 0.
								}
							
							float sum1 = 0;
							
							for (int k=0; k<i; k++)
								{
									sum1 += L[i][k] * U[k][j];
								}
							
							U[i][j] = A[i][j] - sum1;                      // finding the elements above the diagonal of upper triangular matrix from its general expression.
							
						}
					
					for (int i=j+1; i<n; i++)
						{
							U[i][j] = 0;                                   // setting the elements below diagonal of upper triangular matrix as 0.
							
							float sum2 = 0;
							
							for (int k=0; k<j; k++)
								{
									sum2 += L[i][k] * U[k][j];
								}
							
							L[i][j] = (A[i][j] - sum2) / U[j][j];          // finding the elements below diagonal of lower triangular matrix from its general expression. 
						}
				}
			
			float *y;
			y = forward_sub(*L,*B,n);                              // applying forward substitution algorithm to find matrix y where L.y = B. 
			
			float *x;
			x = backward_sub(*U,y,n);                              // applying backward substitution algorithm to find final solution matrix x where U.x = y.
			
			Output = fopen("OutputFile.txt","w");                  // opening the output file.
			
			if (Output == NULL)
				{
					puts("ERROR: In creating OutputFile.txt");
					return 1;
				}
			
			for (int i=0; i<n; i++)
				{
					fprintf(Output,"%f\n",*(x+i));                     // writing the elements of final solution in the output file. 
				}
			
			fclose(Output);
			
			printf("The Solution Matrix is written in OutputFile.txt file\n");
			return 2;
	}
