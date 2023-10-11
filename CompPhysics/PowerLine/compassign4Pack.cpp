#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>

using namespace std;

const int imax = 3;
const int jmax = 3;
const double dx = 1.0 / (imax - 1.0);
const double dy = 1.0 / (jmax - 1.0);

extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

int main() {
    int k;
    int kmax = imax*jmax;
    int kwire = (((imax - 1)/2) * imax + ((jmax -1)/2)); 
    double* RHS = new double[kmax];
    double** A = new double*[kmax];
    double* hold = new double[kmax*kmax];
    int* pivot = new int[kmax];

    int nd, icount, mcount = 0, i, j, i1, j1, ierr, nrhs;

    nd = kmax;
    nrhs = 1;

    for(int ii = 0; ii < kmax; ii++){
        A[ii] = new double[kmax];
    }
    
    for(int ii = 0;ii < kmax;ii++){
        RHS[ii] = 0.0;
    }

    for(int ii = 0;ii < kmax;ii++){
        for(int jj = 0;jj < kmax;jj++){
            A[jj][ii] = 0.0;
        }
    }

    for(int ii = 2; ii < imax;ii++){
        for(int jj = 2; jj < jmax;jj++){
            k = (ii-1) * (imax) + (jj - 1);
            A[k][k] = -4.0;
            A[k+1][k] = 1.0;
            A[k-1][k] = 1.0;
            A[k+imax][k] = 1.0;
            A[k-imax][k] = 1.0;
        }
    }

    for(int ii = 1; ii < imax-1;ii++){
        A[ii][ii] = -4.0;
        A[ii+1][ii] = 1.0;
        A[ii-1][ii] = 1.0;
        A[ii+imax][ii+imax] = 1.0;
        k = (kmax - 1) - ii;
        if(k != kwire){
            A[k][k] = -4.0;
            A[k+1][k] = 1.0;
            A[k-1][k] = 1.0;
            A[k-imax][k] = 1.0;
        }
    }

    for(int jj = 1; jj < jmax-1;jj++){
        k = jmax * jj;
        A[k][k] = -4.0;
        A[k+1][k] = 1.0;
        A[k+imax][k] = 1.0;
        A[k-imax][k] = 1.0;
        k = jmax * jj + (imax-1);
        A[k][k] = -4.0;
        A[k-1][k] = 1.0;
        A[k+imax][k] = 1.0;
        A[k-imax][k] = 1.0;
    }

    A[0][0] = -4.0;
    A[1][0] = 1.0;
    A[imax][0] = 1.0;
    A[imax-1][imax-1] = -4.0;
    A[imax-2][imax-1] = 1.0;
    A[imax+imax-1][imax-imax-1] = 1.0;
    k = (jmax-1) * imax;
    A[k][k] = -4.0;
    A[k+1][k] = 1.0;
    A[k-imax][k] = 1.0;
    k = (jmax-1) * imax + imax - 1;
    A[k][k] = -4.0;
    A[k-1][k] = 1.0;
    A[k-imax][k] = 1.0;

    for(int ii=0;ii < kmax;ii++){
        A[ii][kwire] = 0.0;
    }
    A[kwire][kwire] = 1.0;
    RHS[kwire] = 6.0e4;
    
    dgesv_(&nd,&nrhs,hold,&nd,pivot,RHS,&nd,&ierr);

    for(int i = 0; i < kmax; i++){
        cout << RHS[i] << endl;
    }

    delete[] A, RHS;
}
