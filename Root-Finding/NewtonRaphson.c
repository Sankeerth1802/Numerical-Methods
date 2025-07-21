#include <stdio.h>
#include <math.h>

double function(double x);                                                   		// function to evaluate the given expression.
double f_derivative(double x);                                               		// function to evaluate derivative of given expression.

int main()
	{
		double x_arr[10] = {0};                                                  		// array which stores the value of x.
		
		for (int i=0; i<9; i++)
			{
				x_arr[i+1] = x_arr[i] - function(x_arr[i])/f_derivative(x_arr[i]);  		// Newton Raphson algorithm.
			}
			
		printf("The first positive root of given equation(accurate to 4th decimal) is %.4lf\n",x_arr[9]);
		
		return 0;
	}

double function(double x)
	{
		double x_new = 0;
		x_new = exp(-(x*x)) + 1.5 - (1/(x + 0.2));
		return x_new;
	}

double f_derivative(double x)
	{
		double x_new = 0;
		double square = (x + 0.2) * (x + 0.2);
		x_new = -2 * x * exp(-(x*x)) + 1/square;
		return x_new;
	}
