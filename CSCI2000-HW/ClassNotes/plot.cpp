#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>

using namespace std;

const double PI = 3.14159265359;
double THETA = PI/180; //radian of one degree

int main() {
    const int NUM_SEGMENTS = 360;

    ofstream dataFile;
    dataFile.open("data.dat");

    for (int i = 0; i < NUM_SEGMENTS; i++) {
        double x = i * cos(i * THETA);
        double y = i * sin(i * THETA);

        dataFile << x << " " << y << endl;
    }

    dataFile.close();

    system("gnuplot graph.gnuplot");
    system("display tmp.png");
}


