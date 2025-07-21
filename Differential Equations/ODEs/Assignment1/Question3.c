#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f1(double x, double y1, double y2)
    {
        return cos(x) - y2;
    }

double f2(double x, double y1, double y2)
    {
        return 4*cos(x) - sin(x) - 4*y2 + 3*y1;
    }

void RK4(double *x, double *y1, double *y2, double h, int n)
    {   
        for (int i = 0; i < n - 1; i++)
            {
                double k11 = h * f1(x[i],y1[i],y2[i]);
                double k12 = h * f2(x[i],y1[i],y2[i]);

                double k21 = h * f1(x[i] + h/2, y1[i] + k11/2, y2[i] + k12/2);
                double k22 = h * f2(x[i] + h/2, y1[i] + k11/2, y2[i] + k12/2);

                double k31 = h * f1(x[i] + h/2, y1[i] + k21/2, y2[i] + k22/2);
                double k32 = h * f2(x[i] + h/2, y1[i] + k21/2, y2[i] + k22/2);

                double k41 = h * f1(x[i] + h, y1[i] + k31, y2[i] + k32);
                double k42 = h * f2(x[i] + h, y1[i] + k31, y2[i] + k32);

                y1[i + 1] = y1[i] + (k11 + 2*k21 + 2*k31 + k41)/6;
                y2[i + 1] = y2[i] + (k12 + 2*k22 + 2*k32 + k42)/6;

            }
    }

void ABM4(double *x, double *y1, double *y2, double h, int n)
    {
        RK4(x,y1,y2,h,4);

        for (int i = 3; i < n - 1; i++)
            {
                double y1_p = y1[i] + (h/24)*(55*f1(x[i],y1[i],y2[i]) - 59*f1(x[i - 1],y1[i - 1],y2[i - 1]) + 37*f1(x[i - 2],y1[i - 2],y2[i - 2]) - 9*f1(x[i - 3],y1[i - 3],y2[i - 3]));
                double y2_p = y2[i] + (h/24)*(55*f2(x[i],y1[i],y2[i]) - 59*f2(x[i - 1],y1[i - 1],y2[i - 1]) + 37*f2(x[i - 2],y1[i - 2],y2[i - 2]) - 9*f2(x[i - 3],y1[i - 3],y2[i - 3]));

                y1[i + 1] = y1[i] + (h/24)*(9*f1(x[i + 1],y1_p,y2_p) + 19*f1(x[i],y1[i],y2[i]) - 5*f1(x[i - 1],y1[i - 1],y2[i - 1]) + f1(x[i - 2],y1[i - 2],y2[i - 2]));
                y2[i + 1] = y2[i] + (h/24)*(9*f2(x[i + 1],y1_p,y2_p) + 19*f2(x[i],y1[i],y2[i]) - 5*f2(x[i - 1],y1[i - 1],y2[i - 1]) + f2(x[i - 2],y1[i - 2],y2[i - 2]));
            }
    }

int main() 
    {
        double h = 0.01;  
        double x_final = 5;  
        int n = (int)(x_final / h) + 1; 

        double *x = (double *)malloc(n * sizeof(double));
        double *y1 = (double *)malloc(n * sizeof(double));
        double *y2 = (double *)malloc(n * sizeof(double));

        
        for (int i = 0; i < n; i++) 
            {
                x[i] = i * h;
            }

        y1[0] = 1.0;  
        y2[0] = 0.5;  

        
        ABM4(x, y1, y2, h, n);
        printf("Implementation of ABM4 method is successful.\nThe Solutions obtained are\n");
        
        printf("x\ty1\t\ty2\n");
        for (int i = 0; i < n; i++) 
            {
                printf("%.2f\t%f\t%f\n", x[i], y1[i], y2[i]);
            }

        free(x);
        free(y1);
        free(y2);

        return 0;
    }