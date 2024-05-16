//! @file polyfit.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.
#include <iostream>

void gaussEliminationLS(int m, int n, double** a, double* x) {
    int i, j, k;
    for (i = 0; i < m - 1; i++) {
        //Partial Pivoting
        for (k = i + 1; k < m; k++) {
            //If diagonal element(absolute vallue) is smaller than any of the terms below it
            if (fabs(a[i][i]) < fabs(a[k][i])) {
                //Swap the rows
                for (j = 0; j < n; j++) {
                    double temp;
                    temp = a[i][j];
                    a[i][j] = a[k][j];
                    a[k][j] = temp;
                }
            }
        }
        //Begin Gauss Elimination
        for (k = i + 1; k < m; k++) {
            double  term = a[k][i] / a[i][i];
            for (j = 0; j < n; j++) {
                a[k][j] = a[k][j] - term * a[i][j];
            }
        }

    }
    //Begin Back-substitution
    for (i = m - 1; i >= 0; i--) {
        x[i] = a[i][n - 1];
        for (j = i + 1; j < n - 1; j++) {
            x[i] = x[i] - a[i][j] * x[j];
        }
        x[i] = x[i] / a[i][i];
    }

}

double polyfit(size_t n, size_t deg, const double* xp, const double* yp,
    const double* wp, double* pp)
{
    //arrays to store the c and y-axis data-points
    double* X = new double[2 * deg + 1];
    double* yp_copy = new double[n];
    double* xp_copy = new double[n];
    double del = 1;

    for (int i = 0; i < n; i++) {
        if (wp[0] != -1)
            del = wp[i] * wp[i];
        else
            del = 1;
        yp_copy[i] = yp[i] * del;
        xp_copy[i] = xp[i] / del;
    }

    for (int i = 0; i <= 2 * deg; i++) {
        X[i] = 0;
        for (int j = 0; j < n; j++) {
            if (wp[0] != -1)
                del = wp[j] * wp[j];
            else
                del = 1;
            X[i] = X[i] + pow(xp[j], i) * del;
        }
    }
    //the normal augmented matrix
    double** B = new double* [deg + 1];
    for (int i = 0; i < deg + 1; i++) {
        B[i] = new double[deg + 2];
    }

    // rhs

    double* Y = new double[deg + 1];

    for (int i = 0; i <= deg; i++) {
        Y[i] = 0;
        for (int j = 0; j < n; j++) {
            Y[i] = Y[i] + pow(xp[j], i) * yp_copy[j];
        }
    }

    for (int i = 0; i <= deg; i++) {
        for (int j = 0; j <= deg; j++) {
            B[i][j] = X[i + j];
        }
    }
    for (int i = 0; i <= deg; i++) {
        B[i][deg + 1] = Y[i];
    }
    double* A1 = new double[deg + 1];

    gaussEliminationLS(deg + 1, deg + 2, B, A1);
    for (int i = 0; i <= deg; i++) {
        pp[i] = A1[i];
    }

    double res = 0;
    for (int j = 0; j < n; j++) {
        double yj = 0;
        for (int i = 0; i <= deg; i++)
        {
            yj += pp[i] * pow(xp[j], i);
        }
        res += pow(yj - yp[j], 2);
    }
    delete[] X;
    delete[] yp_copy;
    delete[] xp_copy;
    delete[] A1;
    delete[] Y;
    for (int i = 0; i < deg + 1; i++) {
        delete[] B[i];
    }
    delete[] B;
    //std::cout << "norm fake = " << pow(res / n, 0.5) << "\n\n\n";
    return pow(res / n, 0.5);

}