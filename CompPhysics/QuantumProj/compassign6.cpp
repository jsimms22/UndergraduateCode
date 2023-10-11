#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

const double PI = 4 * 1/tan(1);
const double hBar = 6.62607e-34;
const double L = 4.8e-10; //meters
const int nmax = 100;
const float fnmax = (float) nmax;
const double h = L / (32*nmax);
const double mass = 9.209e-34; //mass of an electron

extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *nrhs, double *w, double *work, int *lwork, int *info);

void WriteData(double** Psis, double** Psic, double* xaxis){
    ofstream output;
    output.open("sin.dat");
    for(int i = 0;i < nmax;i++){
        for(int j = 0;j < nmax;j++){
            output << xaxis[j] << " " << Psis[i][j] << "\n";
        }
    }
    output.close();
    
    output.open("cos.dat");
    for(int i = 0;i < nmax;i++){
        for(int j = 0;j < nmax;j++){
            output << xaxis[j] << " " << Psic[i][j] << "\n";
        }
    }
    output.close();
}

int main(){
    char jobz, uplo;
    double *w, *work;
    int *pivot;
    int numel, numel2, lwork, nrhs, ierr, iterate;
    double sum;
    jobz = 'V'; uplo = 'U';

    numel = nmax;
    numel2 = numel * numel;
    lwork = numel * numel - 1;
    w = new double[numel];
    work = new double[lwork];

    double* Hs = new double[numel2]; double* Hc = new double[numel2];
    double** Psis = new double*[numel]; double** Psic = new double*[numel];
    for(int i = 0;i < numel;i++){
        Psis[i] = new double[numel];
        Psic[i] = new double[numel];
    }
    double* xaxis = new double[nmax];
    for(int i = 0;i < nmax;i++){
        xaxis[i] = -L/2.0 + (L / fnmax) * i; 
    }

    //cout << "seg fault 1?\n";

    for(int n = 1;n <= nmax;n++){
        for(int m = 1;m <= nmax;m++){
            sum = 0.0;
            for(int i = 0;i < 8 *fnmax;i+=4){
                sum += 2.0/45.0 * h * 7.0 * -3.204e-18 * 2.0/L 
                    * sin((n*PI*(1.0+(i/(4.0*fnmax))))/8.0) 
                    * sin((n*PI*(1.0+(i/(4.0*fnmax))))/8.0);
                sum += 2.0/45.0 * h * 32.0 * -3.204e-18 * 2.0/L 
                    * sin((n*PI*(1.0+((i+1.0)/(4.0*fnmax))))/8.0) 
                    * sin((n*PI*(1.0+((i+1.0)/(4.0*fnmax))))/8.0);
                sum += 2.0/45.0 * h * 12.0 * -3.204e-18 * 2.0/L 
                    * sin((n*PI*(1.0+((i+2.0)/(4.0*fnmax))))/8.0) 
                    * sin((n*PI*(1.0+((i+2.0)/(4.0*fnmax))))/8.0);
                sum += 2.0/45.0 * h * 32.0 * -3.204e-18 * 2.0/L 
                    * sin((n*PI*(1.0+((i+3.0)/(4.0*fnmax))))/8.0) 
                    * sin((n*PI*(1.0+((i+3.0)/(4.0*fnmax))))/8.0);
                sum += 2.0/45.0 * h * 7.0 * -3.204e-18 * 2.0/L 
                    * sin((n*PI*(1.0+((i+4.0)/(4.0*fnmax))))/8.0) 
                    * sin((n*PI*(1.0+((i+4.0)/(4.0*fnmax))))/8.0);
            }
            iterate = (n-1)*numel+(m-1);
            Hs[iterate] = 2 * sum;
            iterate = (n-1) + (numel*(m-1));
            Hs[iterate] = 2 * sum;
        }
    }

    //cout << "seg fault 2?\n";

    for(int n = 0;n < nmax;n++){
        for(int m = 0;m < nmax;m++){
            sum = 0.0;
            for(int i = 0;i < 8 *fnmax;i+=4){
                sum += 2.0/45.0 * h * 7.0 * -3.204e-18 * 2.0/L 
                    * cos((n*PI*(1.0+(i/(4.0*fnmax))))/8.0) 
                    * cos((n*PI*(1.0+(i/(4.0*fnmax))))/8.0);
                sum += 2.0/45.0 * h * 32.0 * -3.204e-18 * 2.0/L 
                    * cos((n*PI*(1.0+((i+1.0)/(4.0*fnmax))))/8.0) 
                    * cos((n*PI*(1.0+((i+1.0)/(4.0*fnmax))))/8.0);
                sum += 2.0/45.0 * h * 12.0 * -3.204e-18 * 2.0/L 
                    * cos((n*PI*(1.0+((i+2.0)/(4.0*fnmax))))/8.0) 
                    * cos((n*PI*(1.0+((i+2.0)/(4.0*fnmax))))/8.0);
                sum += 2.0/45.0 * h * 32.0 * -3.204e-18 * 2.0/L 
                    * cos((n*PI*(1.0+((i+3.0)/(4.0*fnmax))))/8.0) 
                    * cos((n*PI*(1.0+((i+3.0)/(4.0*fnmax))))/8.0);
                sum += 2.0/45.0 * h * 7.0 * -3.204e-18 * 2.0/L 
                    * cos((n*PI*(1.0+((i+4.0)/(4.0*fnmax))))/8.0)
                    * cos((n*PI*(1.0+((i+4.0)/(4.0*fnmax))))/8.0);
            }
            iterate = n * numel + m;
            Hc[iterate] = 2 * sum;
            iterate = n + numel * m;
            Hc[iterate] = 2 * sum;
        }
    }

    //cout << "seg fault 3?\n";

/*    for(int i = 0;i < nmax;i++){
        for(int j = 0;j < nmax;j++){
            cout << setprecision(4) << Hs[j*numel+i] << " ";
        }
        cout << "\n";
    }*/

    for(int i = 0;i < nmax;i++){
        int n = i + 1;
        iterate = i * numel + i;
        Hs[iterate] += (hBar * hBar * n * n * PI * PI)/(2 * mass * L * L);
    }

    for(int n = 0;n < nmax;n++){
        iterate = n * numel + n;
        Hc[iterate] += (hBar * hBar * n * n * PI * PI)/(2 * mass * L * L);
    }

/*    for(int i = 0;i < nmax;i++){
        for(int j = 0;j < nmax;j++){
            cout << setprecision(4) << Hs[j*numel+i] << " ";
        }
        cout << "\n";
    }*/

    //cout << "seg fault 4?\n";

    dsyev_(&jobz,&uplo,&numel,Hs,&numel,w,work,&lwork,&ierr);

    dsyev_(&jobz,&uplo,&numel,Hc,&numel,w,work,&lwork,&ierr);


/*    for(int i = 0;i < nmax;i++){
        for(int j = 0;j < nmax;j++){
            cout << setprecision(4) << Hs[j*numel+i] << " ";
        }
        cout << "\n";
    }*/

    cout << ierr << "\n";
    double x = 0.0;
    for(int nx = 0;nx < nmax;nx++){
        x = nx * h - (L/2);
        for(int m = 1;m <= nmax;m++){
            for(int n = 1;n <= nmax;n++){
                iterate = (m-1)*numel+(n-1);
                Psis[m-1][nx] = Hs[iterate] * sin((n * PI * x)/L);
            }
        }
    }
    x = 0.0;
    for(int nx = 0;nx < nmax;nx++){
        x = nx * h - (L/2.0);
        for(int m = 0;m < nmax;m++){
            for(int n = 0;n < nmax;n++){
                iterate = m*numel+n;
                Psic[m][nx] = Hc[iterate] * cos((n * PI * x)/L);
            }
        }
    }
/*    cout << "\n";
    cout << "\n";
    cout << "Psi below:\n";
    for(int i = 0;i < nmax;i++){
        for(int j = 0;j < nmax;j++){
            cout << setprecision(4) << Psis[i][j] << " ";
        }
        cout << "\n";
    }*/

    //cout << "seg fault 5?\n";


    WriteData(Psis,Psic,xaxis);

    //cout << "seg fault 6?\n";
    
    delete[] Hs;
    delete[] Hc;
    delete[] xaxis;
    delete[] Psis;
    delete[] Psic;
    delete w;
    delete work;
    delete pivot;
}










