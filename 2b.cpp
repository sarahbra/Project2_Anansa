#include<iostream>
#include <cmath>
#include "A.h"

void jacobi_method (double ** A, double ** R, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 0; j++) {
            if (i == 0) {
                R[i][j] = 1.0;
            } else {
                R[i][j] = 0.0;
            }
        }
    }

    int row, col;
    double epsilon = 1.0*10^(-8);
    int iterations = 0;
    int max_iterations = n*n*n;
    double max_nondiag = maxnondiag(A, &col, &row, n);

    while (iterations < max_iterations && fabs(max_nondiag) > epsilon) {
        max:nondiag = maxnondiag(A, &col, &row, n);
        rotate (A, r, col, row, n);
        iterations++;
    }
    cout << "Number of iterations: " << iterations << "\n";
    return;
}

void rotate (double ** A, double ** R, int col, int row, int n) {
    double c, s;
    if(A[row][col] != 0) {

    }
}

void maxnondiag (double ** A, int * col, int * row, int n) {
    double max = 0.0;
    int temp;
    for (int i=0; i<n; i++) {
        for (int j=i+1; j<n; j++){
            temp = fabs(A[i][j]);
            if (temp > max) {
                max = temp;
                *row = i;
                *col = j;
            }
        }
    }
}


