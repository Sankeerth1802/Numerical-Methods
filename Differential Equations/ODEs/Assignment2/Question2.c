#include <stdio.h>
#include <stdlib.h>
#include "matrix_inverse.h"

#define tolerance 1e-6

void find_F(double *y, double *F, double h, int n)
    {
        // function that contains all the discretised form of equations.

        y[0] = 1.0;                                 // boundary values
        y[n + 1] = 0.0;

        for (int i = 1; i < n + 1; i++)
            {
                // discretising at all the internal points.

                double y_xx = (y[i + 1] - 2*y[i] + y[i - 1])/(h*h);             // second derivative
                double y_x = (y[i + 1] - y[i - 1]) / (2*h);                     // first derivative

                F[i - 1] = (1 + y[i])*y_xx + (y_x * y_x);
            }
        
    }

void evaluate_jacobian(double *y, double *J, double h, int n)
    {
        // function to find the jacobian matrix

        for (int i = 0; i < n; i++)
            {
                // intitialising all the elements to 0

                for (int j = 0; j < n; j++)
                    {
                        J[i*n + j] = 0.0;
                    }
            }
        
        J[0*(n) + 0] = (y[2] - 4*y[1] + y[0] - 2)/(h*h);                   // values of the first row
        J[0*(n) + 1] = (2 + 2*y[1] + y[2] - y[0])/(2*h*h);

        for (int i = 2; i < n; i++)
            {
                // for all the rows between 0 and n

                J[(i - 1)*(n) + (i - 2)] = (2 + 2*y[i] - y[i + 1] + y[i - 1])/(2*h*h);
                J[(i - 1)*(n) + (i - 1)] = (y[i + 1] - 4*y[i] + y[i - 1] - 2)/(h*h);
                J[(i - 1)*(n) + i] = (2 + 2*y[i] + y[i + 1] - y[i - 1])/(2*h*h);

            }
        
        J[(n - 1)*(n) + (n - 2)] = (2 + 2*y[n] - y[n + 1] + y[n - 1])/(2*h*h);         // last row values
        J[(n - 1)*(n) + (n - 1)] = (y[n + 1] - 4*y[n] + y[n - 1] - 2)/(h*h);
        
    }

int main()
    {
        double h = 0.01;                                                      
        int n = (int)(1/h) - 1;                                               // number of internal points

        double *y = (double *)malloc((n + 2)*sizeof(double));                 // y values including the boundary values.
        double *F = (double *)malloc(n*sizeof(double));                       // F functions only for internal points.
        double *J = (double *)malloc((n*n)*sizeof(double));                   // Jacobian only at internal points
        double *delta = (double *)malloc(n * sizeof(double));                 // to store j_inv times F
        
        double *J_inv = (double *)calloc(n*n,sizeof(double)); 

        y[0] = 1.0;                                           // boundary values
        y[n + 1] = 0.0;

        for (int i = 1; i < (n + 1); i++)
            {
                // initial guesses

                y[i] = 1.0;
            }
        
        int iter_count = 0;
        int max_iter = 50;

        while (iter_count < max_iter)
            {
                // implementing multi variate newton raphson method

                find_F(y,F,h,n);
                evaluate_jacobian(y,J,h,n);
                J_inv = invert_matrix(J,n,n);

                for (int i = 0; i < n; i++)
                    {
                        // finding the delta vector

                        delta[i] = 0.0;

                        for (int j = 0; j < n; j++)
                            {
                                delta[i] -= J_inv[i*n + j] * F[j];
                            }
                    }

                for (int i = 1; i < (n + 1); i++)
                    {
                        // updating y values

                        y[i] += delta[i - 1];
                    }

                iter_count++;
                
                double F_norm = 0.0;

                for (int i = 0; i < n; i++) 
                    {
                        // calculating the 1 norm(sum of absolute values)

                        F_norm += fabs(F[i]);
                    }

                if (F_norm < tolerance) 
                    {
                        // convergence checking.

                        printf("Solutions are found using the Newton Raphson method.\n");
                        break;
                    }
            }
        
        if (iter_count == max_iter) 
            {
                // divergence case

                printf("Failed to converge.\n");
            } 

        printf("The solutions at every interval of 10 nodes are:\nx\ty\n");

        for (int i = 0; i <= n+1; i += 10) 
            {
                printf("%.2f\t%.6f\n", i*h, y[i]);
            }
        
        free(J_inv);
        free(y);
        free(F);
        free(J);
        free(delta);

        return 0;
    }