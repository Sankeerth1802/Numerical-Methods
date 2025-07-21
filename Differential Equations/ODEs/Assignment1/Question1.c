#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f1(double u1, double u2)
    {
        return u2;
    }

double f2(double u1, double u2)
    {
        return -25*u2 - u1;
    }

void RK2(double *t, double *u1, double *u2, double h, int n)
    {
        for (int i = 0; i < n - 1; i++)
            {
                double k11 = h * f1(u1[i],u2[i]);
                double k12 = h * f2(u1[i],u2[i]);

                double k21 = h * f1(u1[i] + k11,u2[i] + k12);
                double k22 = h * f2(u1[i] + k11,u2[i] + k12);

                u1[i + 1] = u1[i] + (k11 + k21)/2;
                u2[i + 1] = u2[i] + (k12 + k22)/2;
            }
    }

void ABM2(double *t, double *u1, double *u2, double h, int n)
    {
        RK2(t,u1,u2,h,2);

        for (int i = 1; i < n - 1; i++)
            {
                double u1_p = u1[i] + (h/2)*(3*f1(u1[i],u2[i]) - f1(u1[i - 1],u2[i - 1]));
                double u2_p = u2[i] + (h/2)*(3*f2(u1[i],u2[i]) - f2(u1[i - 1],u2[i - 1]));

                u1[i + 1] = u1[i] + (h/2)*(f1(u1[i],u2[i]) + f1(u1_p,u2_p));
                u2[i + 1] = u2[i] + (h/2)*(f2(u1[i],u2[i]) + f2(u1_p,u2_p));
            }
    }

int main()
    {
        double h[2] = {0.1, 0.12};
        int x_final = 1;

        for (int i = 0; i < 2; i++)
            {
                int n = (int)(x_final/h[i]) + 1;

                double *x = (double *)malloc(n * sizeof(double));
                double *u1 = (double *)malloc(n * sizeof(double));
                double *u2 = (double *)malloc(n * sizeof(double));
                double *u3 = (double *)malloc(n * sizeof(double));
                double *u4 = (double *)malloc(n * sizeof(double));

                for (int j = 0; j < n; j++)
                    {
                        x[j] = j * h[i];
                    }
                
                u1[0] = 1;
                u2[0] = 0;
                u3[0] = 1;
                u4[0] = 0;

                printf("\nThe Solution obtained using RK2 for h = %.2f are:\n",h[i]);

                RK2(x,u1,u2,h[i],n);

                printf("x\tu1\t\tu2\n");

                for (int j = 0; j < n; j++)
                    {
                        printf("%.2f\t%f\t%f\n",x[j],u1[j],u2[j]);
                    }
                
                printf("\nThe Solution obtained using ABM2 for h = %.2f are:\n",h[i]);

                ABM2(x,u3,u4,h[i],n);

                printf("x\tu1\t\tu2\n");

                for (int j = 0; j < n; j++)
                    {
                        printf("%.2f\t%f\t%f\n",x[j],u3[j],u4[j]);
                    }
                
                free(x);
                free(u1);
                free(u2);
                free(u3);
                free(u4);
            }
        
        return 0;
    }