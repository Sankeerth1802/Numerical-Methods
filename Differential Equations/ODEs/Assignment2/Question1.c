#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define tolerance 1e-6

double f0(double u0, double u1)
    {
        // first derivative of y.
        return u1;
    }

double f1(double u0, double u1)
    {
        // second derivative of y.
        return -(u1 * u1)/(1 + u0);
    }

double error(double u0,double u_actual)
    {
        // to calculate the error of the boundary value at the end.
        return (u_actual - u0);
    }

void RK2(double *x, double *u0, double *u1, double h, int n)
    {
        // 2nd order rk function

        for (int i = 0; i < n - 1; i++)
            {
                double k11 = h * f0(u0[i],u1[i]);
                double k12 = h * f1(u0[i],u1[i]);

                double k21 = h * f0(u0[i] + k11, u1[i] + k12);
                double k22 = h * f1(u0[i] + k11, u1[i] + k12);

                u0[i + 1] = u0[i] + (k11 + k21)/2;
                u1[i + 1] = u1[i] + (k12 + k22)/2;
            }
    }

int main()
    {
        double h = 0.01;                                           // step size
        int n = (int)(1/h) + 1;                                    // total number of points.

        double *x = (double *)malloc(n * sizeof(double));
        double *u0 = (double *)malloc(n * sizeof(double));
        double *u1 = (double *)malloc(n * sizeof(double));

        for (int i = 0; i < n; i++)
            {
                // initialising the x values

                x[i] = i * h;
            }
        
        u0[0] = 1;                                                 // initial value of y
        double u0_actual = 0;                                      // outer boundary value of y

        double u1_0 = -0.1;                                        // 1st guess of initial slope 
        double u1_1 = 0.1;                                         // 2nd guess
        double u1_updated;                                         // updated value of initial slope 

        int maximum_iterations = 100;                              // defining maximum number of iterations and iteration count.
        int iteration_count = 0;

        while (iteration_count != maximum_iterations)
            {
                // implementing secant method to find the initial slope

                u1[0] = u1_0;
                RK2(x,u0,u1,h,n);
                double err1 = error(u0[n - 1],u0_actual);          // calculating the error using 1st guess

                u1[0] = u1_1;
                RK2(x,u0,u1,h,n);
                double err2 = error(u0[n - 1],u0_actual);          // calculating error with 2nd guess

                u1_updated = u1_1 - (err2*(u1_1 - u1_0))/(err2 - err1);  // updating the initial slope using the guesses and errors.

                if (fabs(err2) < tolerance)
                    {
                        // convergence condition

                        printf("The value of u'(0) is successfully found using the Secant Method.\n");
                        break;
                    }
                
                u1_0 = u1_1;                                   // updating the values.
                u1_1 = u1_updated;
                iteration_count++;
            }
        
        if (iteration_count == maximum_iterations)
            {
                // divergence case

                printf("ERROR: Appropriate value for u'(0) so that it satisfies u(1) = 0 is not found.\n");
            }

        else
            {
                u1[0] = u1_1;
                RK2(x,u0,u1,h,n);                                           // finally implementing rk with the found initial slope.

                printf("The values of x,y and y' at an interval of 10 are:\n%-12s %-12s %-12s\n", "x", "y", "y'");

                for (int i = 0; i < n; i += n/10) 
                    {
                        printf("%-12.2f %-12f %-12f\n", x[i], u0[i], u1[i]);
                    }
            }
        
        free(x);
        free(u0);
        free(u1);

        return 0;
    }