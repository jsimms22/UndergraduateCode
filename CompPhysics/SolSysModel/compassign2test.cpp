#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

const double tau = 60 * 60 * 24;
const double G = 6.67408e-11;
const int nbodies = 10;
const int ntime = 256 * 365;

double Accel(int n, double** x, double** y, double** z, double m[nbodies], double radius, double Ax[nbodies], double Ay[nbodies], double Az[nbodies]) {
    for(int i = 0; i < nbodies; i++) {
        Ax[i] = 0;
        Ay[i] = 0;
        Az[i] = 0;
        for(int j = 0; j < nbodies; j++) {
            if (i != j) {
                radius = sqrt(pow(x[j][n] - x[i][n],2) + pow(y[j][n] - y[i][n],2) 
                    + pow(z[j][n] - z[i][n],2));
                Ax[i] = ((G * m[j] * (x[j][n] - x[i][n])) / pow(radius, 3))
                    + Ax[i];
                Ay[i] = ((G * m[j] * (y[j][n] - y[i][n])) / pow(radius, 3)) 
                    + Ay[i];
                Az[i] = ((G * m[j] * (z[j][n] - z[i][n])) / pow(radius, 3)) 
                    + Az[i];
            }
        }
    }
    
    return Ax[nbodies], Ay[nbodies], Az[nbodies];
}

double Momentum() {

}

void writeData(double** x, double** y, double** z){

    string filename1 = "data.dat";
    ofstream dataFile1;
    dataFile1.open(filename1.c_str());

    for(int b = 0; b < nbodies; b++){
        for(int n = 0; n < ntime - 1; n++) {
            dataFile1 << x[b][n] << " " << y[b][n] << " " << z[b][n] << endl;;
        }
    }
        
    dataFile1.close();

    for(int b = 0; b < nbodies; b++){
        stringstream ss;
        ss << b;
        string str;
        ss >> str;
        string filename = "data" + str + ".dat";

        ofstream dataFile;
        dataFile.open(filename.c_str());
        for(int n = 0; n < ntime - 1; n++) {
            dataFile << x[b][n] << " " << y[b][n] << " " << z[b][n] << endl;
        }
        dataFile.close();
    }
    string write;
    for(int b = 0; b < nbodies; b++){
        stringstream ss;
        ss << b;
        string str;
        ss >> str;
        if (b != 0) {
            write = write + ", \'data" + str + ".dat\'";
        } else {
            write = "\'data" + str + ".dat\'";
        }
    }  

    /*string plot = "graph.gnuplot";

    ofstream gnuIchooseYou;
    gnuIchooseYou.open(plot.c_str());

    gnuIchooseYou << "set terminal png medium size 900,600 background '#FFFFFF'" << endl;
    gnuIchooseYou << "set xlabel \"x\"" << endl;
    gnuIchooseYou << "set ylabel \"y\"" << endl;
    gnuIchooseYou << "set key on" << endl;
    gnuIchooseYou << "set output \"tmp.png\"" << endl;
    gnuIchooseYou << "plot \"data0.dat\" with lines, \"data1.dat\" with lines" , \"data2.dat\" with lines, \"data3.dat\" with lines, \"data4.dat\" with lines ", \"data5.dat\" with lines, \"data6.dat\" with lines, \"data7.dat\" with lines, \"data8.dat\" with lines, \"data9.dat\"" << endl;

    gnuIchooseYou.close();
    string pngFilename = "mercury.png";

    system("gnuplot graph.gnuplot");
    rename("tmp.png", pngFilename.c_str());
    string write2 = "display " + pngFilename;
    system(write2.c_str());*/

}

int main() {
    double **x = new double*[nbodies];
    double **y = new double*[nbodies];
    double **z = new double*[nbodies];
    double **Vx = new double*[nbodies]; 
    double **Vy = new double*[nbodies];
    double **Vz = new double*[nbodies];

    for(int i = 0; i < nbodies; i++) {
        x[i] = new double[ntime];
        y[i] = new double[ntime];
        z[i] = new double[ntime];
        Vx[i] = new double[ntime];
        Vy[i] = new double[ntime];
        Vz[i] = new double[ntime];
    }

    //Sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune, pluto
    double m[nbodies] = {1988550000e21, 330.2e21, 4868.5e21, 5973.6e21, 641.85e21, 1898600e21, 568460e21, 86832e21, 102430e21, 13.105e21};
    double radius = 0;
    double Ax[nbodies], Ay[nbodies], Az[nbodies];
    //double Lx[nbodies][ntime], Ly[nbodies][ntime], Lz[nbodies][ntime];
    //double LxSum[ntime], LySum[ntime], LzSum[nbodies];

    //Sun
    x[0][0] = 5.388e5 * 1000; y[0][0] = 4.906e5 * 1000; z[0][0] = -2.402e4 * 1000;
    Vx[0][0] = -3.16e-3 * 1000; Vy[0][0] = 1.192e-2 * 1000; Vz[0][0] = 6.211e-5 * 1000;
    //Mercury
    x[1][0] = 2.320e7 * 1000; y[1][0] = 4.088e7 * 1000; z[1][0] = 1.198e6 * 1000;
    Vx[1][0] = -5.221e1 * 1000; Vy[1][0] = 2.578e1 * 1000; Vz[1][0] = 6.895e0 * 1000;
    //Venus
    x[2][0] = -2.638e7 * 1000; y[2][0] = -1.047e8 * 1000; z[2][0] = 8.679e4 * 1000;
    Vx[2][0] = 3.369e1 * 1000; Vy[2][0] = -8.813e0 * 1000; Vz[2][0] = -2.065e0 * 1000;
    //Earth
    x[3][0] = 1.503e8 * 1000; y[3][0] = 8.615e6 * 1000; z[3][0] = -2.479e4 * 1000;
    Vx[3][0] = -2.093 * 1000; Vy[3][0] = 2.965e1 * 1000; Vz[3][0] = -2.041e-3 * 1000;
    //Mars
    x[4][0] = 1.473e8 * 1000; y[4][0] = -1.465e8 * 1000; z[4][0] = -6.707e6 * 1000;
    Vx[4][0] = 1.806e1 * 1000; Vy[4][0] = 1.921e1 * 1000; Vz[4][0] = -4.090e-2 * 1000;
    //Jupiter
    x[5][0] = -8.134e8 * 1000; y[5][0] = -4.747e7 * 1000; z[5][0] = 1.839e7 * 1000;
    Vx[5][0] = 6.097e-1 * 1000; Vy[5][0] = -1.243e1 * 1000; Vz[5][0] = 3.760e-2 * 1000;
    //Saturn
    x[6][0] = -3.537e8 * 1000; y[6][0] = -1.459e9 * 1000; z[6][0] = 3.944e7 * 1000;
    Vx[6][0] = 8.857 * 1000; Vy[6][0] = -2.306 * 1000; Vz[6][0] = -3.125e-1 * 1000;
    //Uranus
    x[7][0] = 2.766e9 * 1000; y[7][0] = 1.121e9 * 1000; z[7][0] = -3.167e7 * 1000;
    Vx[7][0] = -2.708e0 * 1000; Vy[7][0] = 5.994e0 * 1000; Vz[7][0] = 5.587e-2 * 1000;
    //Neptune
    x[8][0] = 4.225e9 * 1000; y[8][0] = -1.493e9 * 1000; z[8][0] =  -6.663e7 * 1000;
    Vx[8][0] = 1.774e0 * 1000; Vy[8][0] = 5.157e0 * 1000; Vz[8][0] = -1.468e-1 * 1000;
    //Pluto
    x[9][0] = 1.401e9 * 1000; y[9][0] = -4.761e9 * 1000; z[9][0] = 1.042e8 * 1000;
    Vx[9][0] = 5.343e0 * 1000; Vy[9][0] = 4.162e-1 * 1000; Vz[9][0] = -1.606e0 * 1000;
    
    for (int b = 0; b < nbodies; b++) {
        Accel(0, x, y, z, m, radius, Ax, Ay, Az);
        Vx[b][1] = Vx[b][0] + Ax[b] * tau;
        Vy[b][1] = Vy[b][0] + Ay[b] * tau;
        Vz[b][1] = Vz[b][0] + Az[b] * tau;
        x[b][1] = Vx[b][1] * tau + x[b][0];
        y[b][1] = Vy[b][1] * tau + y[b][0];
        z[b][1] = Vz[b][1] * tau + z[b][0];
    }
    
    int j, n;
    for(n = 1; n < ntime - 1; n++) {
        Accel(n, x, y, z, m, radius, Ax, Ay, Az);
        for(j = 0; j < nbodies; j++) {
            x[j][n+1] = 2 * x[j][n] - x[j][n-1] + Ax[j] * tau * tau;
            y[j][n+1] = 2 * y[j][n] - y[j][n-1] + Ay[j] * tau * tau;
            z[j][n+1] = 2 * z[j][n] - z[j][n-1] + Az[j] * tau * tau;
            Vx[j][n+1] = (x[j][n+1] - x[j][n]) / tau;
            Vy[j][n+1] = (y[j][n+1] - y[j][n]) / tau;
            Vz[j][n+1] = (z[j][n+1] - z[j][n]) / tau;
        }
    }

    /*for(int n = 0; n < ntime - 1; n++){
        cout << "Step: " << n << endl;
        cout << x[8][n] << " " << y[8][n] << " " << z[8][n] << endl;
    }*/

    writeData(x, y, z);
    delete x, y, z, Vx, Vy, Vz;
}
