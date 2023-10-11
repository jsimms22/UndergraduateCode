#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>

using namespace std;

const int nx = 101;
const int ny = 101;
const double dx = (1.0 / (nx - 1.0));
const double dy = (1.0 / (ny - 1.0));
const int nsteps = 4001;
const double tau = (1.0e-2 / (nsteps - 1.0));
const double R = 0.25;
const double rho = 1;
const double velocity = 343;
const double c = velocity / rho;

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

double WriteData(double*** z){
    for(int t = 0; t < nsteps; t++){
        if( t % 50 == 0){
            stringstream ss;
            ss << t;
            string str;
            ss >> str;
            string filename = "time" + str + ".dat";
            ofstream output;
            output.open(filename.c_str());
    
            for(int i = 0; i < nx; i++){
                for(int j = 0; j < ny; j++){
                    output << i*dx << " " << j*dy << " "<< z[i][j][t] << endl;
                }
            }

            output.close();
        }
    }
}

int main() {
    double cx = ((nx - 1) / 2);
    double cy = ((ny - 1) / 2);

    //cout << 1 << endl;

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

    //cout << 2 << endl;

    //************* dynamically allocating a - acceleration 2dimensional array
    double **a = new double*[nx];
    for(int i = 0; i < nx; i++){
        a[i] = new double[ny];
    }
    //*************

    //cout << 3 << endl;

    //Initializing z to be 0 at all points in the field at all times observed
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            for(int t = 0; t < nsteps; t++){
                z[i][j][t] = 0;
                vz[i][j][t] = 0;
            }
        }
    }

    //cout << 4 << endl;
    
    //Initializing a to be 0 at all points in the field at a desired t - time
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            a[i][j] = 0;
        }
    }
    
    //cout << 5 << endl;

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

    //cout << 6 << endl;
    
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

    //cout << 7 << endl;

    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            a1[i][j] = 0; a2[i][j] = 0; a3[i][j] = 0; a4[i][j] = 0;
            k1[i][j] = 0; k2[i][j] = 0; k3[i][j] = 0; k4[i][j] = 0;

        }
    }

    //cout << 8 << endl;

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

        //cout << "t = " << t << " ended" << endl;
    }

    //cout << 8 << endl;

    WriteData(z);

    delete[] z, vz, a, a1, a2, a3, a4, k1, k2, k3, k4;

    return 0;
}
