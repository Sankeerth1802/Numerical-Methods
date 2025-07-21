#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define GAUSS_SIEDEL_H

double* GaussSiedel(double *C, double *B,long int x,double tol)
	{ 
		double r[x][1], rm, *y;
		
		y = (double*)malloc(x*sizeof(double));
		
		double y_prev[x][1];
		
		for (int i=0; i<x; i++)
			{
				*(y+i) = 1;
			}
		do
			{
				for (int i = 0; i < x; i++)
							{
								y_prev[i][0] = *(y+i);
							}
							
				for (int i = 0; i < x; i++)
					{
							
						double sum = 0;
						
						for (int j = 0; j < x; j++)
							{
								if (j != i)
									{
										sum += *(C+(i*x)+j) * *(y+j);
									}
							}
							
						*(y+i) = (*(B+i) - sum) / *(C + i*(x+1));
					
					}
						rm = 0;
						
						for (int i = 0; i < x; i++)
							{
								r[i][0] = *(y+i) - y_prev[i][0];
								
								rm += r[i][0] * r[i][0];
							}
							
						rm = pow(rm, 0.5);
						
					
					
			} while(rm > tol);
			
		return y;
	}

double* ConjugateGradient(double *C,double *B,long int x,double tol)
	{
		double y_prev[x][1], *y, r_0[x][1], r[x][1], d_0[x][1], d[x][1], rmax;
		
		y = (double*)malloc(x*sizeof(double));
		
		for (int i = 0; i < x; i++)
			{
				*(y+i) = 1;
			}
			
		double mid1[x][1];
		
		for (int i=0; i<x; i++)
			{
				mid1[i][0] = 0;
				
				for (int k=0; k<x; k++)
					{
						mid1[i][0] += *(C+(i*x)+k) * *(y+k);  
					}
			}
			
		for (int i=0; i<x; i++)
			{
				r_0[i][0] = *(B+i) - mid1[i][0];
			}	
				
		for (int i=0; i<x; i++)
			{
				d_0[i][0] = r_0[i][0];
				
				r[i][0] = r_0[i][0];
			}
		do
			{	
				for (int i = 0; i < x; i++)
					{
						y_prev[i][0] = *(y+i);
					}
					
				double num=0,den=0;
				
				for (int i=0; i<x; i++)
					{
						num += r_0[i][0] * d_0[i][0];
					}
					
				double mid[x][1];
				
				for (int i=0; i<x; i++)
					{	
						mid[i][0] = 0;
						
						for(int k=0; k<x; k++)
							{
								mid[i][0] += *(C+(i*x)+k) * d_0[k][0];
							}
					}
					
				for (int i=0; i<x; i++)
					{
						den += d_0[i][0] * mid[i][0];
					}
				
				double alpha = num/den;
				
				for (int i=0; i<x; i++)
					{
						*(y+i) = y_prev[i][0] + (alpha*d_0[i][0]);
						
						r[i][0] = r_0[i][0] - (alpha*mid[i][0]);
					}
					
				double num1=0;
				
				for (int i=0; i<x; i++)
					{
						num1 += r[i][0] * mid[i][0];
					}
					
				double beta = -(num1)/den;
				
				for (int i=0; i<x; i++)
					{
						d[i][0] = r[i][0] + (beta*d_0[i][0]);
						
						r_0[i][0] = r[i][0];
						
						d_0[i][0] = d[i][0];
					}
					
				rmax = 0;
				
				for (int i=0; i<x; i++)
					{ 
						rmax += r_0[i][0] * r_0[i][0];
					}
					
				rmax = pow(rmax,0.5); 
				
			} while (rmax > tol);
		return y;
	}
