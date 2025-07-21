#include <stdio.h>
#include <stdlib.h>
#define M_PI 3.14159265358979323846
#include <math.h>


double k = 0.02;
double m = 1;
double g = 9.807;

double f1(double u1, double u2, double u3, double u4)
    {
        return u3;
    }

double f2(double u1, double u2, double u3, double u4)
    {
        return u4;
    }

double f3(double u1, double u2, double u3, double u4)
    {
        return -(k/m)*u3*sqrt(pow(u3,2) + pow(u4,2));
    }

double f4(double u1, double u2, double u3, double u4)
    {
        return  -g - (k/m)*u4*sqrt(pow(u3,2) + pow(u4,2));
    }

void RK4(double *t, double *u1, double *u2, double *u3, double *u4, double h, int n)
    {
        for (int i = 0; i < n - 1; i++)
            {
                double k11 = h * f1(u1[i],u2[i],u3[i],u4[i]);
                double k12 = h * f2(u1[i],u2[i],u3[i],u4[i]);
                double k13 = h * f3(u1[i],u2[i],u3[i],u4[i]);
                double k14 = h * f4(u1[i],u2[i],u3[i],u4[i]);

                double k21 = h * f1(u1[i] + k11/2,u2[i] + k12/2,u3[i] + k13/2,u4[i] + k14/2);
                double k22 = h * f2(u1[i] + k11/2,u2[i] + k12/2,u3[i] + k13/2,u4[i] + k14/2);
                double k23 = h * f3(u1[i] + k11/2,u2[i] + k12/2,u3[i] + k13/2,u4[i] + k14/2);
                double k24 = h * f4(u1[i] + k11/2,u2[i] + k12/2,u3[i] + k13/2,u4[i] + k14/2);

                double k31 = h * f1(u1[i] + k21/2,u2[i] + k22/2,u3[i] + k23/2,u4[i] + k24/2);
                double k32 = h * f2(u1[i] + k21/2,u2[i] + k22/2,u3[i] + k23/2,u4[i] + k24/2);
                double k33 = h * f3(u1[i] + k21/2,u2[i] + k22/2,u3[i] + k23/2,u4[i] + k24/2);
                double k34 = h * f4(u1[i] + k21/2,u2[i] + k22/2,u3[i] + k23/2,u4[i] + k24/2);

                double k41 = h * f1(u1[i] + k31,u2[i] + k32,u3[i] + k33,u4[i] + k34);
                double k42 = h * f2(u1[i] + k31,u2[i] + k32,u3[i] + k33,u4[i] + k34);
                double k43 = h * f3(u1[i] + k31,u2[i] + k32,u3[i] + k33,u4[i] + k34);
                double k44 = h * f4(u1[i] + k31,u2[i] + k32,u3[i] + k33,u4[i] + k34);

                u1[i + 1] = u1[i] + (k11 + 2*k21 + 2*k31 + k41)/6;
                u2[i + 1] = u2[i] + (k12 + 2*k22 + 2*k32 + k42)/6;
                u3[i + 1] = u3[i] + (k13 + 2*k23 + 2*k33 + k43)/6;
                u4[i + 1] = u4[i] + (k14 + 2*k24 + 2*k34 + k44)/6;

            }
        
    }

int main() 
    {
        double v0 = 20.0;  
        double theta = 30.0 * M_PI / 180.0; 
        double h = 0.01;  
        int n = 10000;  

        double *t = (double *)malloc(n * sizeof(double));
        double *u1 = (double *)malloc(n * sizeof(double));
        double *u2 = (double *)malloc(n * sizeof(double));
        double *u3 = (double *)malloc(n * sizeof(double));
        double *u4 = (double *)malloc(n * sizeof(double));

        
        t[0] = 0.0;
        u1[0] = 0.0;  
        u2[0] = 0.0;  
        u3[0] = v0 * cos(theta);  
        u4[0] = v0 * sin(theta);  

        RK4(t, u1, u2, u3, u4, h, n);
        printf("Implementation of RK4 is successful.\n");

        int i;

        for (i = 1; i < n; i++) 
            {
                t[i] = t[i-1] + h;

                if (u2[i] < 0) 
                    {  
                        break;
                    }
            }

        double time = t[i-1]; 
        double distance = u1[i-1];

        printf("Time of flight: %.4f seconds\n", t[i - 1]);
        printf("Distance travelled along x before striking the ground is: %f meters\n", u1[i - 1]);

        free(t);
        free(u1);
        free(u2);
        free(u3);
        free(u4);

        return 0;
    }