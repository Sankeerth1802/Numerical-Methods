#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x, double t)
    {
        return (4*x)/t + pow(t,4) * exp(t);
    }

double *Solution(double *time, double h, int n)
    {
        double *y = (double *)malloc(n * sizeof(double));

        for (int i = 0; i < n; i++)
            {
                *(y + i) = pow(*(time + i),4) * (exp(*(time + i)) - exp(1));
            }
        
        return y;
    }

double *RK4(double *time, double h, int n)
    {   
        double *y = (double *)malloc(n * sizeof(double));

        for (int i = 0; i < n; i++)
            {
                *(y + i) = 0;
            }
        
        for (int i = 0; i < n - 1; i++)
            {
                double k1 = h * f(y[i],time[i]);
                double k2 = h * f(y[i] + k1/2,time[i] + h/2);
                double k3 = h * f(y[i] + k2/2,time[i] + h/2);
                double k4 = h * f(y[i] + k3,time[i] + h);

                y[i + 1] = y[i] + (k1 + 2*k2 + 2*k3 + k4)/6;
            }
        
        return y;

    }

double *ABM4(double *y, double *time, double h, int n)
    {
        double *y1 = (double *)malloc(n * sizeof(double));
        
        for (int i = 0; i < 4; i++)
            {
                y1[i] = y[i];
            }
        
        for (int i = 3; i < n - 1; i++)
            {
                double y_p = y1[i] + (h/24)*(55*f(y1[i],time[i]) - 59*f(y1[i - 1],time[i - 1]) + 37*f(y1[i - 2],time[i - 2]) - 9*f(y1[i - 3],time[i - 3]));

                y1[i + 1] = y1[i] + (h/24)*(9*f(y_p,time[i + 1]) + 19*f(y1[i],time[i]) - 5*f(y1[i - 1],time[i - 1]) + f(y1[i - 2],time[i - 2])); 
            }
        
        return y1;
    }

int main()
    {
        double h = 0.1;
        int n = (int)(1/h) + 1;

        double *time = (double *)malloc(n * sizeof(double));
        int *err = (int *)malloc(n * sizeof(int));
        int *err1 = (int *)malloc(n * sizeof(int));
        
        time[0] = 1;
        
        for (int j = 1; j < n; j++)
            {
                time[j] = time[j - 1] + h;
            }

        double *sol = RK4(time, h, n);
        double *truesol = Solution(time, h, n);
        double *sol1 = ABM4(sol, time, h, n);

        printf("The solutions obtained for h = %.2f are given below\n",h);
        printf("%-10s %-15s %-15s %-15s %-15s %-15s\n", "Time", "True Solution", "RK4 Solution", "ABM4 Solution", "RK4 Errors", "ABM4 Errors");

        for (int i = 0; i < n; i++)
            {
                err[i] = fabs(sol[i] - truesol[i])*1e6;
                err1[i] = fabs(sol1[i] - truesol[i])*1e6;

                printf("%-10.2f %-15.6f %-15.6f %-15.6f %-15d %-15d\n",time[i],truesol[i],sol[i],sol1[i],err[i],err1[i]);                
            }

        printf("\n");

        free(time);
        free(sol);
        free(truesol);
        free(sol1);
        free(err);
        free(err1);

        return 0;
    }
