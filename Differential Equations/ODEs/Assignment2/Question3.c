#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_inverse.h"

#define PI 3.14159265358979323846
#define tolerance 1e-6

void find_F(double *y, double *F, double h, int n)
    {
        // function that contains the discretised forms od given equation
        
        F[0] = (2*(y[1] - y[0]))/(h*h) - 1;                // for the first node, using ghost node; neumann bc gave u_{-1} = u_1

        for (int i = 1; i < n + 1; i++)
            {
                // for all the nodes in between the boundary values

                double r = 1.0 + i*h;                               // current value of r
                double u_r = (y[i+1] - y[i-1])/(2*h);               // first derivative
                double u_rr = (y[i+1] - 2*y[i] + y[i-1])/(h*h);     // the second derivative
                
                F[i] = u_rr + u_r/r - cos(PI*(r-1.0)/2.0);
            }

        double u_ghost = y[n] + 0.2*h*(pow(y[n + 1],4) - 0.4);      // for the last node, applying ghost; robin bc gave relation between u_n+2 and u_n and u_n+1
        double u_rr = (u_ghost - 2*y[n + 1] + y[n])/(h*h);
        double u_r = (u_ghost - y[n])/(2*h);
        double r = 1.0 + (n + 1)*h;

        F[n + 1] = u_rr + u_r/r;
    }

void evaluate_jacobian(double *y, double *J, double h, int n)
    {
        // calculating the jacobian matrix

        for (int i = 0; i < (n + 2); i++)
            {
                // initialising all the elements to 0

                for (int j = 0; j < (n + 2); j++)
                    {
                        J[i*(n + 2) + j] = 0.0;
                    }
            }

        J[0*(n + 2) + 0] = -2/(h*h);                             // the first row of the jacobian matrix
        J[0*(n + 2) + 1] = 2/(h*h);

        for (int i = 1; i < (n + 1); i++)
            {
                // all the rows in between 0(r = 1.0) and n+2(r = 2.0)

                double r = 1.0 + i*h;

                J[i*(n + 2) + (i - 1)] = 1/(h*h) - 1/(2*h*r);
                J[i*(n + 2) + i] = -2/(h*h);
                J[i*(n + 2) + (i + 1)] = 1/(h*h) + 1/(2*h*r);
            }
        
        double r = 1.0 + (n + 1)*h;                             // the last row of jacobian
        J[(n + 1)*(n + 2) + n] = 2/(h*h);
        J[(n + 1)*(n + 2) + (n + 1)] = (0.8*h*pow(y[n + 1],3) - 2)/(h*h) + (0.4*pow(y[n + 1],3))/r;
    }

int main()
    {
        double h = 0.1;  
        int n = (int)((2.0-1.0)/h) - 1;                           // here again, n represents the number of internal points
        
        
        double *y = (double *)malloc((n + 2)*sizeof(double));     // size of all the matrices is same as no boundary value is known. 
        double *F = (double *)malloc((n + 2)*sizeof(double));
        double *J = (double *)malloc(((n + 2)*(n + 2))*sizeof(double));
        double *delta = (double *)malloc((n + 2)*sizeof(double));
        double *J_inv = (double *)malloc(((n + 2)*(n + 2))*sizeof(double));
        
    
        for (int i = 0; i < (n + 2); i++)
            {
                // initial guess for y values

                y[i] = -1.0;
            }
        
        int iter_count = 0;
        int max_iter = 50; 
        
        while (iter_count < max_iter)
            {   
                
                find_F(y, F, h, n);
                evaluate_jacobian(y, J, h, n);
                J_inv = invert_matrix(J, n + 2, n + 2);

                if (J_inv == NULL)
                    {
                        // error handling

                        printf("Matrix inversion failed!\n");
                        break;
                    }
                
                for (int i = 0; i < n + 2; i++)
                    {
                        // finding the product of jacobian inverse and F

                        delta[i] = 0.0;
                        for (int j = 0; j < n + 2; j++)
                            {
                                delta[i] -= J_inv[i*(n + 2) + j] * F[j];
                            }
                    }
                

                for (int i = 0; i < (n + 2); i++)
                    {
                        // updating the values of y

                        y[i] += delta[i];
                    }
                
                double F_norm = 0.0;

                for (int i = 0; i < n + 2; i++)
                    {
                        F_norm += fabs(F[i]);
                    }
                
                if (F_norm < tolerance)
                    {
                        // convergence check

                        printf("Solutions are successfully obtained using Newton Raphson method.\n");
                        break;
                    }
                
                iter_count++;
            }
        
        printf("Solutions for every interval of 10 nodes are:\nr\tu(r)\n");

        for (int i = 0; i < (n + 2); i += (n + 2)/10)
            {
                double r = 1.0 + i*h;
                printf("%.2f\t%f\n", r, y[i]);
            }
        
        free(y);
        free(F);
        free(J);
        free(delta);
        free(J_inv);
        
        return 0;
    }