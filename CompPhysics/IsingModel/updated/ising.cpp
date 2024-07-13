#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>

// new headers
#include <memory>
#include <matplot/matplot.h>
#include "matrix.hpp"

constexpr int N = 50;

int main(){
    std::ofstream output;
    output.open("isingPara.dat");

    srand(time(NULL));
    auto mu = std::make_shared<matrix::Matrix2<double,N,N>>();
    
    double E, delE, ESqr, ESum, EAvg, ESqrAvg, mSum; // E = energy
    long int m, mSqr, bigNum; // m = total observable electrons
    double mAvg, mSqrAvg;
    double J, muB, B, k, kT; // see below for explanations
    int temp1, temp2; double temp3; // random value holders
    int tmp1, tmp2; // random lattice position holders
    int i1, i2, j1, j2;
    double mVar, EVar;

    J = 1.0; // angular momentum quantum number
    B = 20.0; // magnetic field strength, Tesla
    muB = 5.788381e-5; // bohr mangneton, eV per Tesla
    k = 8.6173303e-5; // boltzmann's constant, eV per Kelvin
    m = N*N; // total observable electrons in lattice
    bigNum = 100000;

    std::vector<double> X{};
    std::vector<double> Y{};

    for(int i = 0;i < N;i++){
        for(int j = 0;j < N;j++){
            (*mu)[i][j] = 1; // (*mu) will always be +/- 1 at every point
        }
    }
    int i, j;

    for(int iT = 0; iT < 2500; iT++){
        kT = k * iT;
        E = 0.0;
        mSum = 0;

        for(i = 0;i < N;i++){
            for(j = 0;j < N;j++){
                E += -B * muB * (*mu)[i][j];
            }
        }

        for(i = 0;i < N;i++){
            for(j = 0;j < N;j++){
                i1 = i - 1; i2 = i + 1;
                j1 = i - 1; j2 = j + 1;
                if(i1 >= 0){
                    E += -J * (*mu)[i][j] * (*mu)[i1][j];
                } if(i2 < N){
                    E += -J * (*mu)[i][j] * (*mu)[i2][j];
                } if(j1 >= 0){
                    E += -J * (*mu)[i][j] * (*mu)[i][j1];
                } if(j2 < N){
                    E += -J * (*mu)[i][j] * (*mu)[i][j2];
                }
            }
        }

        for(int val = 0;val < bigNum;val++){
            temp1 = rand()%50;
            temp2 = rand()%50;
            
            tmp1 = temp1; tmp2 = temp2;
            (*mu)[tmp1][tmp2] = -1 * (*mu)[tmp1][tmp2];

            delE = -B * muB * (*mu)[tmp1][tmp2] * 2;
            i1 = tmp1 - 1; i2 = tmp2 + 1;
            j1 = tmp2 - 1; j2 = tmp2 + 1;
            i = tmp1; j = tmp2;
            
            if(i1 >= 0){
                delE += -J * (*mu)[i][j] * (*mu)[i1][j] * 2;
            } if(i2 < N){
                delE += -J * (*mu)[i][j] * (*mu)[i2][j] * 2;
            } if(j1 >= 0){
                delE += -J * (*mu)[i][j] * (*mu)[i][j1] * 2;
            } if(j2 < N){
                delE += -J * (*mu)[i][j] * (*mu)[i][j2] * 2;
            }
            //cout << delE << endl;

            temp3 = (rand()%1000000) / 1.0e6;
            if((delE > 0) && (temp3 > exp(-delE / kT))){
                (*mu)[tmp1][tmp2] = -1 * (*mu)[tmp1][tmp2];
                //m += 2*(*mu)[tmp1][tmp2];
                //E += delE;

            }else{
                m += 2*(*mu)[tmp1][tmp2];
                E += delE;
                //(*mu)[tmp1][tmp2] = -1 * (*mu)[tmp1][tmp2];

            }

            ESqr += E * E;
            mSqr += m * m;
            ESum += E;
            mSum += m;
            //cout << ESqr << " " << mSqr << endl;
        }

        EAvg = ESum / bigNum;
        mAvg = (mSum / bigNum);

        X.push_back(iT);
        Y.push_back(mAvg);

        output << iT << " " << mAvg << std::endl;

        // ESqrAvg = ESqr / bigNum;
        // mSqrAvg = (double)(mSqr / bigNum);

        // EVar = ESqrAvg - EAvg * EAvg;
        // mVar = mSqrAvg - mAvg * mAvg;

        // std::cout << kT << ", " << EAvg << ", " << sqrt(EVar) << ", " << mAvg << ", " << sqrt(mVar) << "\n";
        
        std::cout << iT << '\n';
    }

    matplot::scatter(X,Y);
    matplot::show();
    output.close();
}
