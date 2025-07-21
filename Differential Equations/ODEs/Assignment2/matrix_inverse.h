#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *invert_matrix(double *matrix, int row, int col) 
    {
        double *aug = (double *)malloc(row * 2 * col * sizeof(double));         // augumented matrix
        
        for (int i = 0; i < row; i++) 
            {
                for (int j = 0; j < 2 * col; j++) 
                    {
                        if (j < col)
                            {
                                aug[i * (2 * col) + j] = matrix[i * col + j];   // first part of augumented matrix.
                            }
                            
                        else
                            {
                                if (i == (j - col))
                                    {
                                        aug[i * (2 * col) + j] = 1.0;           // initialising the second part as identity matrix
                                    }
                                
                                else
                                    {
                                        aug[i * (2 * col) + j] = 0.0;
                                    }
                            }
                    }
            }
        
        for (int i = 0; i < row; i++) 
            {
                // converting first part to identity matrix by implementing pivoted gaussian elimination.

                double pivot = aug[i * (2 * col) + i];                     // pivot element.
                if (fabs(pivot) < 1e-10) 
                    {
                        // case when diagonal elements are 0.

                        free(aug);
                        return NULL;
                    }
                
                for (int j = 0; j < 2 * col; j++) 
                    {
                        // scaling the entire row by the pivot element

                        aug[i * (2 * col) + j] /= pivot;
                    }
                
                for (int k = 0; k < row; k++) 
                    {
                        // eliminating remaining elements in the present column.

                        if (k != i) 
                            {
                                double factor = aug[k * (2 * col) + i];     // multiplication coefficient

                                for (int j = 0; j < 2 * col; j++) 
                                    {
                                        aug[k * (2 * col) + j] -= factor * aug[i * (2 * col) + j];
                                    }
                            }
                    }
            }
        
        double *inverse = (double *)malloc(row * col * sizeof(double));      // inverse which is the second part of augumented matrix.

        for (int i = 0; i < row; i++) 
            {
                for (int j = 0; j < col; j++) 
                    {
                        inverse[i * col + j] = aug[i * (2 * col) + (j + col)];   // calculating inverse separately from the augumented matrix
                    }
            }
        
        free(aug);

        return inverse;
    }
