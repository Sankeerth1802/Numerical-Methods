#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define STEPHEN_CONST 5.670374419e-8
#define tolerance 1e-6

double rho = 7900;
double d = 0.002;
double T_f = 1500;

double f(double T)
    {
        double k = (2 * STEPHEN_CONST)/(rho * d);
        return (k*(pow(T_f,4) - pow(T,4))) / (0.324*T + 446.47);
    }

void RK4(double *t, double *T, double h, int n)
    {
        for (int i = 0; i < n - 1; i++)
            {
                double k1 = h * f(T[i]);
                double k2 = h * f(T[i] + k1/2);
                double k3 = h * f(T[i] + k2/2);
                double k4 = h * f(T[i] + k3);

                T[i + 1] = T[i] + (k1 + 2*k2 + 2*k3 + k4)/6;
            }
    }

int main()
    {
        double h = 0.01; 
        int n = 100000; 

        double *t = malloc(n * sizeof(double));
        double *T = malloc(n * sizeof(double));

        t[0] = 0;
        T[0] = 300; 

        for (int j = 0; j < n; j++)
            {
                t[j] = j * h;
            }

        RK4(t, T, h, n);
        printf("Implementation of RK4 is successful.\n");
        
        int i;
        for (i = 0; i < n; i++)
        {
            if (fabs(T[i] - T_f) < tolerance)
            {
                printf("Time taken to reach thermal equilibrium: %.2f seconds\n", t[i]);
                break;
            }
        }

        free(t);
        free(T);
        return 0;
    }