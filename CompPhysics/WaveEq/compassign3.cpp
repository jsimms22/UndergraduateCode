#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>

using namespace std;

constexpr int nx = 101;
constexpr int ny = 101;
constexpr double dx = (1.0 / (nx - 1.0));
constexpr double dy = (1.0 / (ny - 1.0));
constexpr int nsteps = 4001;
constexpr double tau = (1.0e-2 / (nsteps - 1.0));
constexpr double R = 0.25;
constexpr double rho = 1;
constexpr double velocity = 343;
constexpr double c = velocity / rho;

double Acceleration(int t, double*** z, double** a){
    double dx2, dy2;
    for(int i = 1; i < nx - 1; i++){
        for(int j = 1; j < ny -1; j++){
            dx2 = (z[i+1][j][t] - 2 * z[i][j][t] + z[i-1][j][t]) / (dx * dx);
            dy2 = (z[i][j+1][t] - 2 * z[i][j][t] + z[i][j-1][t]) / (dy * dy);
            a[i][j] = c * c * (dx2 + dy2);
        }
    }
    return** a;
}

double WaveEquation(int i,int j,int t){
    return exp(-((i * dx) * (i * dx) + (j * dy) * (j * dy)) / (R * R)) - exp(-1);
}

void WriteData(double*** z){
    for(int t = 0; t < nsteps; t+=50){
        // std::cout << t << "\n";
        if( t % 50 == 0){
            stringstream ss;
            ss << t;
            string str;
            ss >> str;
            string filename = "time" + str + ".dat";
            // had to replace ofstream with fopen, ofstream was crashing due to size
            FILE* output = fopen(filename.c_str(),"w+");
    
            for(int i = 0; i < nx; i++){
                std::string line{};
                for(int j = 0; j < ny; j++){
                    line = std::to_string(i*dx) + " " + std::to_string(j*dy) + " " + std::to_string(z[i][j][t]) + "\n";
                    fputs(line.c_str(),output);
                }
            }
            fclose(output);
        }
    }
}

int main() {
    double cx = ((nx - 1) / 2);
    double cy = ((ny - 1) / 2);

    //************* dynamically allocating z - position 3dimensional array
    double ***z = new double**[nx];
    for(int i = 0; i < nx; i++){ 
        z[i] = new double*[ny];
        for(int j = 0; j < ny; j++){ 
            z[i][j] = new double[nsteps];
        }
    }
    //*************
    
    double ***vz = new double**[nx];
    for(int i = 0; i < nx; i++){
        vz[i] = new double*[ny];
        for(int j = 0; j < ny; j++){
            vz[i][j] = new double[nsteps];
        }
    }

    //************* dynamically allocating a - acceleration 2dimensional array
    double **a = new double*[nx];
    for(int i = 0; i < nx; i++){
        a[i] = new double[ny];
    }
    //*************

    //Initializing z to be 0 at all points in the field at all times observed
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            for(int t = 0; t < nsteps; t++){
                z[i][j][t] = 0;
                vz[i][j][t] = 0;
            }
        }
    }

    //Initializing a to be 0 at all points in the field at a desired t - time
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            a[i][j] = 0;
        }
    }

    double r;
    int i1, i2, j1, j2;
    for(int i = 0; i < cx / 2; i++){
        for(int j = 0; j < cy / 2; j++){
            i1 = i + cx; i2 = -i + cx;
            j1 = j + cy; j2 = -j + cy;
            r = sqrt((i * dx) * (i * dx) + (j * dy) * (j * dy));
            if(R >= r) {
                z[i1][j1][0] = WaveEquation(i,j,0);
                z[i1][j2][0] = WaveEquation(i,j,0);
                z[i2][j1][0] = WaveEquation(i,j,0);
                z[i2][j2][0] = WaveEquation(i,j,0);
            } else {
                z[i1][j1][0] = 0;
                z[i1][j2][0] = 0;
                z[i2][j1][0] = 0;
                z[i2][j2][0] = 0;
            }
        }
    }
    
    double **a1 = new double*[nx]; double **a2 = new double*[nx];
    double **a3 = new double*[nx]; double **a4 = new double*[nx];

    double **k1 = new double*[nx]; double **k2 = new double*[nx]; 
    double **k3 = new double*[nx]; double **k4 = new double*[nx];

    for(int i = 0; i < nx; i++){
        a1[i] = new double[ny]; a2[i] = new double[ny]; a3[i] = new double[ny];
            a4[i] = new double[ny];

        k1[i] = new double[ny]; k2[i] = new double[ny]; k3[i] = new double[ny]; 
            k4[i] = new double[ny];
    }

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            a1[i][j] = 0; a2[i][j] = 0; a3[i][j] = 0; a4[i][j] = 0;
            k1[i][j] = 0; k2[i][j] = 0; k3[i][j] = 0; k4[i][j] = 0;

        }
    }

    for(int t = 0; t < nsteps; t++){
        //cout << "t = " << t << endl;
        Acceleration(t,z,a);
        for(int i = 0; i < nx - 1; i++){
            for(int j = 0; j < ny - 1; j++){
                vz[i][j][t+1] = vz[i][j][t] + a[i][j] * tau;
                k1[i][j] = vz[i][j][t] + a[i][j] * tau;
                z[i][j][t+1] = z[i][j][t] + k1[i][j] * tau;
            }
        }
        Acceleration(t,z,a1);
        for(int i = 1; i < nx - 1; i++){
            for(int j = 1; j < ny - 1; j++){
                
                k2[i][j] = vz[i][j][t+1] + a1[i][j] * (tau / 2.0);
                z[i][j][t+1] = z[i][j][t] + k2[i][j] * (tau / 2.0);
            }
        }
        Acceleration(t,z,a2);
        for(int i = 1; i < nx - 1; i++){
            for(int j = 1; j < ny - 1; j++){
                k3[i][j] = vz[i][j][t+1] + a2[i][j] * (tau / 2.0);
                z[i][j][t+1] = z[i][j][t] + k3[i][j] * (tau / 2.0);
            }
        }
        Acceleration(t,z,a3);
        for(int i = 1; i < nx - 1; i++){
            for(int j = 1; j < ny - 1; j++){
                k4[i][j] = vz[i][j][t+1] + a3[i][j] * tau;
                z[i][j][t+1] = z[i][j][t] + k4[i][j] * tau;
            }
        }
        Acceleration(t,z,a4);
        for(int i = 1; i < nx - 1; i++){
            for(int j = 1; j < ny -1; j++){
                vz[i][j][t+1] = vz[i][j][t] + 1.0/6.0 * (a1[i][j] + 2 * a2[i][j] + 
                    2 * a3[i][j] + a4[i][j]) * tau;
                z[i][j][t+1] =  z[i][j][t] + 1.0/6.0 * (k1[i][j] + 2 * k2[i][j] +
                    2 * k3[i][j] + k4[i][j]) * tau;
            }
        }

    }
    WriteData(z);

    delete[] z, vz, a, a1, a2, a3, a4, k1, k2, k3, k4;

    return 0;
}
