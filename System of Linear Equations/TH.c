#include <stdio.h>
#include "gauss_siedel.h"

int main()
	{
		FILE *ptr,*ptr1;                                                            // to write temperatures
	  int m,n,Ts,Tb;
	  printf("Enter the number of rows\n");
	  scanf("%d",&m);
	  printf("Enter the number of columns\n");
	  scanf("%d",&n);
	  printf("Enter the temperature of bottom points\n");
	  scanf("%d",&Tb);
	  printf("Enter the temperature of sides and top points\n");
	  scanf("%d",&Ts);
	  
	  long int x;
	  x = (m - 1) * (n - 1);
	  
	  double C[x][x] ;                                      // Coefficient matrix
	  double T[x][1] ;                                      // Temperature (Variable) column matrix
	  double B[x][1] ;                                      // Constant matrix
	   
	  for (int i=0; i<x; i++)
	  	{
	  		for (int j=0; j<x; j++)
	  			{
	  				C[i][j] = 0;
	  			}
	  	} 
	  
	  for (int i=0; i<x; i++)
	  	{
	  		T[i][0] = 0;
	  		B[i][0] = 0;	
	  	} 	
		
	
	   
	  for (int j=1; j<n; j++)                                             // Finding Coefficients of Coefficient matrix by converting internal points as 1D array.
	  	{ 
	  		for (int i=1; i<m; i++)                                         // (i,j)th point is expressed as kth point where k = (j - 1) * (m - 1) + i. i runs from 1 to m-1 and j runs from 1 to n-1.
	  			{
	  				if (j == 1)
	  					{
	  						if (i == 1)
	  							{
	  								C[0][0] = 4;
	  								C[0][1] = -1;
	  								C[0][m-1] = -1;
	  								B[0][0] = Ts+Tb;
	  							}
	  						
	  					  else if (i == m-1)
	  							{
	  								C[m-2][m-2] = 4;
	  								C[m-2][m-3] = -1;
	  								C[m-2][2 * m - 3] = -1;
	  								B[m-2][0] = Ts+Tb;
	  							}
	  						
	  						else 
	  							{
	  								C[i-1][i-1] = 4;
	  								C[i-1][i] = -1;
	  								C[i-1][i-2] = -1;
	  								C[i-1][i+m-2] = -1;
	  								B[i-1][0] = Tb;
	  							}	
	  					}
	  				
	  				else if (j == n-1)
	  					{
	  						if (i == 1)
	  							{
	  								C[x-m+1][x-m+1] = 4;
	  								C[x-m+1][x-m+2] = -1;
	  								C[x-m+1][x-(2*m)+2] = -1;
	  								B[x-m+1][0] = Ts+Ts;
	  							}
	  						
	  						else if (i == m-1)
	  							{
	  								C[x-1][x-1] = 4;
	  								C[x-1][x-2] = -1;
	  								C[x-1][x-m] = -1;
	  								B[x-1][0] = Ts+Ts;
	  							}
	  						
	  						else
	  							{
	  								C[x-m+i][x-m+i] = 4;
	  								C[x-m+i][x-m+i+1] = -1;
	  								C[x-m+i][x-m+i-1] = -1;
	  								C[x-m+i][x-(2*m)+i+1] = -1;
	  								B[x-m+i][0] = Ts;
	  							}
	  					}
	  					
	  				else if (i == 1 && j != 1 && j != n-1)
	  					{
	  						long int z = 0;
	  						z = (j-1)*(m-1) + 1;
	  						C[z-1][z-1] = 4;
	  						C[z-1][z] = -1;
	  						C[z-1][z+m-2] = -1;
	  						C[z-1][z-m] = -1;
	  						B[z-1][0] = Ts;
	  					}
	  				
	  				else if (i == m-1 && j != 1 && j != n-1)
	  					{
	  						long int z=0;
	  						z = (j)*(m-1);
	  						C[z-1][z-1] = 4;
	  						C[z-1][z-2] = -1;
	  						C[z-1][z+m-2] = -1;
	  						C[z-1][z-m] = -1;
	  						B[z-1][0] = Ts;
	  					}
	  				
	  				else
	  					{
	  						long int z = 0;
	  						z = (j-1)*(m-1) + i;
	  						C[z-1][z-1] = 4;
	  						C[z-1][z] = -1;
	  						C[z-1][z-2] = -1;
	  						C[z-1][z+m-2] = -1;
	  						C[z-1][z-m] = -1;
	  					}
	  			}
	  	}	
	
		double *Sol;
		Sol = GaussSiedel(*C,*B,x,0.000001);                           // Gauss Siedel algorithm implementation.
		
		ptr = fopen("GaussSiedel.txt","w");
		
		for (int j=n-1; j>0; j--)
			{
				int k = 0;
				
				for (int i=1; i<m; i++)
					{
						k = (j - 1) * (m - 1) + i;
						
						fprintf(ptr,"%lf ",*(Sol+k-1));
					}
				fprintf(ptr,"\n");
			}
		
		fclose(ptr);
		printf("Temperatures of internal points obtained by Gauss Siedel are written in GaussSiedel.txt\n");
		
		double *Sol1;
		Sol1 = ConjugateGradient(*C,*B,x,0.000001);                    // Conjugate Gradient algorithm implementation.
		
		ptr1 = fopen("ConjugateGradient.txt","w");
		
		for (int j=n-1; j>0; j--)
			{
				int k = 0;
				
				for (int i=1; i<m; i++)
					{
						k = (j - 1) * (m - 1) + i;
						
						fprintf(ptr1,"%lf ",*(Sol1+k-1));
					}
				fprintf(ptr1,"\n");
			}
		
		fclose(ptr1);
		printf("Temperatures of internal points obtained by Conjugate Gradient are written in ConjugateGradient.txt\n");
		return 0;
	}


