#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <omp.h>

using namespace std;

const int imax = 25;
const int jmax = 25;
const int kmax = imax * jmax;
const double dx = 1.0 / (imax - 1.0);
const double dy = 1.0 / (jmax - 1.0);

/*void WriteData(double** physical){
    ofstream output;
    output.open("voltmap.dat");

    for(int i = 0; i < imax+2; i++){
        for(int j = 0; j < jmax+2; j++){
            output << physical[i][j] << " ";
        }
        output << endl;
    }

    output.close();   
}*/

int main() {
    clock_t t;
    t = clock();

    int k;
    double ierror;
    int kwire, kwirecenter; 
    double big = 0.0;
    double* temp = new double[kmax];
    int index = 0;
    double factor;
    double tempRHS;
    double* RHS = new double[kmax];
    double** A = new double*[kmax];
    for(int i = 0; i < kmax; i++){
        A[i] = new double[kmax];
    }

    bool** logic = new bool*[imax+2];
    for(int j = 0;j < jmax+2; j++){
        logic[j] = new bool[jmax+2];
    }

    double** physical = new double*[imax+2];
    for(int j = 0;j < jmax+2; j++){
        physical[j] = new double[jmax+2];
    }

    for(int i = 0; i < imax+2;i++){
        for(int j = 0; j < jmax+2; j++){
            physical[i][j] = 0.0;
        }
    }
    for(int i = 1; i < imax+1;i++){
        for(int j = 1; j < jmax+1;j++){
            physical[i][j] = 1.0;
        }
    }
    physical[(imax+1)/2][(jmax+1)/2] = 2.0;

    for(int i = 0;i < imax+2;i++){
        for(int j = 0;j < jmax+2;j++){
            if(physical[i][j] == 1){
                logic[i][j] = false;
            } else {
                logic[i][j] = true;
            }
        }
    }

    for(int i = 0;i < kmax;i++){
        RHS[i] = 0.0;
    }

    for(int i = 0;i < kmax;i++){
        for(int j = 0;j < kmax;j++){
            A[i][j] = 0.0;
        }
    }
//void omp_set_thread_limit(2);
//#pragma omp parallel for schedule(runtime)
for(int z = 0;z < 60;z++){
    int kwire = (((imax - 1)/2) * imax + ((jmax -1)/2));
    //int kwiremoving = dx * dx * z * z + (((imax - 1)/4) * imax + ((jmax-1)/4));

    cout << z << "\n";
    for(int i = 2; i < imax;i++){
        for(int j = 2; j < jmax;j++){
            k = (i-1) * (imax) + (j - 1);
            A[k][k] = -4.0;
            A[k][k+1] = 1.0;
            A[k][k-1] = 1.0;
            A[k][k+imax] = 1.0;
            A[k][k-imax] = 1.0;
        }
    }

    for(int i = 1; i < imax-1;i++){
        A[i][i] = -4.0;
        A[i][i+1] = 1.0;
        A[i][i-1] = 1.0;
        A[i][i+imax] = 1.0;
        k = (kmax - 1) - i;
        if(k != kwire){
            A[k][k] = -4.0;
            A[k][k+1] = 1.0;
            A[k][k-1] = 1.0;
            A[k][k-imax] = 1.0;
        }
    }

    for(int j = 1; j < jmax-1;j++){
        k = jmax * j;
        A[k][k] = -4.0;
        A[k][k+1] = 1.0;
        A[k][k+imax] = 1.0;
        A[k][k-imax] = 1.0;
        k = jmax * j + (imax-1);
        A[k][k] = -4.0;
        A[k][k-1] = 1.0;
        A[k][k+imax] = 1.0;
        A[k][k-imax] = 1.0;
    }

    A[0][0] = -4.0;
    A[0][1] = 1.0;
    A[0][imax] = 1.0;
    A[imax-1][imax-1] = -4.0;
    A[imax-1][imax-2] = 1.0;
    A[imax-1][imax+imax-1] = 1.0;
    k = (jmax-1) * imax;
    A[k][k] = -4.0;
    A[k][k+1] = 1.0;
    A[k][k-imax] = 1.0;
    k = (jmax-1) * imax + imax - 1;
    A[k][k] = -4.0;
    A[k][k-1] = 1.0;
    A[k][k-imax] = 1.0;
    for(int i=0;i < kmax;i++){
        A[kwire][i] = 0.0;
    }
    A[kwire][kwire] = 1.0;
    RHS[kwire] = 6.0e4;

//BOTTOM OF MATRIX CANCELING
    for(int n = 0;n < kmax;n++){
    #pragma omp parallel for private(factor) schedule(runtime)
        for(int j = n+1; j < kmax; j++){
            factor = A[j][n] / A[n][n];
            for(int i = 0; i < kmax; i++){
                A[j][i] = A[j][i] - factor*A[n][i];
            }
            RHS[j] = RHS[j] - factor*RHS[n];
        }
   #pragma end parallel for
    }

//Dividing each equation by the diagonal number to achieve RREF
    double tempnum = 0;
    for(int i = 0; i < kmax; i++){
        tempnum = A[i][i];
        for(int j = 0; j < kmax; j++){
            A[i][j] = A[i][j] / tempnum;
        }
        RHS[i] = RHS[i] / tempnum;
    }

    RHS[kmax-1] = RHS[kmax-1] / A[kmax-1][kmax-1];
    for(int n = kmax-2;n > -1; n--){
        for(int i = n+1; i <  kmax; i++){
            RHS[n] = RHS[n] - RHS[i] * A[n][i];
        }
        RHS[n] = RHS[n] / A[n][n];
    }
}
//#pragma end parallel for
    t = clock() - t;

    cout << ((float)t)/CLOCKS_PER_SEC << endl;

    k = 0;
    for(int i = 1; i < imax+1; i++){
        for(int j = 1; j < jmax+1; j++){
            physical[i][j] = RHS[k];
            k++;
        }
    }

    //WriteData(physical);

    delete[] A, RHS, temp, physical, logic;
}
