//Name: Jeremiah Simmons
//Description: Program to generate plot of different types
//Assignment: CSCI 2000 Assignment 5

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>

using namespace std;

const double PI = 3.14149265359;
const double THETA = PI / 90; //radian step of two degrees

void getCurve(int &curveNum);

void getPngFilename(string &pngFilename);

void getNumPoints(int &numPoints);

void circle(double xValues[], double yValues[], int numPoints, double r);

void epicycloid(double xValues[], double yValues[], int numPoints, double R,
double r);

void hypotrochoid(double xValues[], double yValues[], int numPoints, 
double R, double r, double d);

void hypocycloid(double xValues[], double yValues[], int numPoints,
double R, double r);

void butterfly(double xValues[], double yValues[], int numPoints);

void writeData(double xValues[], double yValues[], int numPoints,
string pngFilename);

int main() {
    
    int curveNum;
    getCurve(curveNum);
    //cout << curveNum << endl;

    string pngFilename;
    getPngFilename(pngFilename);
    //cout << pngFilename << endl;

    int numPoints;
    getNumPoints(numPoints);
    //cout << numPoints << endl;

    double xValues[numPoints], yValues[numPoints], r, R, d;

    switch (curveNum) {
        case 1:
            cout << "Enter data for the circle" << endl;
            cout << "Enter r: ";
            cin >> r;
            circle(xValues, yValues, numPoints, r);

            /*for (int i = 0;i < numPoints; i++) {
                cout << xValues[i] << " " << yValues[i] << endl;
            }*/
        break;
        case 2:
            cout << "Enter data for the epicycloid" << endl;
            cout << "Enter r: ";
            cin >> r;
            cout << "Enter R: ";
            cin >> R;
            epicycloid(xValues, yValues, numPoints, R, r);

           /* for (int i = 0;i < numPoints; i++) {
                cout << xValues[i] << " " << yValues[i] << endl;
            }*/
        break;
        case 3:
            cout << "Enter data for the hypotrochoid" << endl;
            cout << "Enter r: ";
            cin >> r;
            cout << "Enter R: ";
            cin >> R;
            cout << "Enter d: ";
            cin >> d;
            hypotrochoid(xValues, yValues, numPoints, R, r, d);

            /*for (int i = 0; i < numPoints; i++) {
                cout << xValues[i] << " " << yValues[i] << endl;
            }*/
        break;
        case 4:
            cout << "Enter date for the hypocycloid" << endl;
            cout << "Enter r: ";
            cin >> r;
            cout << "Enter R: ";
            cin >> R;
            hypocycloid(xValues, yValues, numPoints, R, r);
            
            /*for (int i = 0; i < numPoints; i++) {
                cout << xValues[i] << " " << yValues[i] << endl;
            }*/
        break;
        case 5:
            butterfly(xValues, yValues, numPoints);
        break;
    }
    
    writeData(xValues, yValues, numPoints, pngFilename);

}

void getCurve(int &curveNum){
    cout << "Which curve to draw?\n"
         << "   1) Circle\n"
         << "   2) Epicycloid\n"
         << "   3) Hypotrochoid\n"
         << "   4) Hypocycloid\n"
         << "   5) Butterfly\n"
         << "Choice: ";    
    cin >> curveNum;
}

void getPngFilename(string &pngFilename){
    cout << "Image name to write to: ";
    cin >> pngFilename;
}

void getNumPoints(int &numPoints){
    cout << "How many points to generate: ";
    cin >> numPoints;
    
    while (numPoints < 2) {
        cout << "Enter integer value of 2 or higher" << endl;
        cout << "How many points to generate: ";
        cin >> numPoints;
    }
}

void circle(double xValues[], double yValues[], int numPoints, double r){
    for (int i = 0; i < numPoints; i++) {
        xValues[i] = r * cos(i * THETA);
        yValues[i] = r * sin(i * THETA);
    }
}

void epicycloid(double xValues[], double yValues[], int numPoints, double R,
double r){
    for (int i = 0; i < numPoints; i++) {
        xValues[i] = (R + r) * cos(i * THETA) - (r * cos(((R + r) / r) *
        (i * THETA)));
        yValues[i] = (R + r) * sin(i * THETA) - (r * sin(((R + r) / r) * 
        (i * THETA)));
    }
}

void hypotrochoid(double xValues[], double yValues[], int numPoints, 
double R, double r, double d){
    for (int i = 0; i < numPoints; i++) {
        xValues[i] = (R - r) * cos(i * THETA) + d * cos(((R - r) / r) * 
        (i * THETA));
        yValues[i] = (R - r) * sin(i * THETA) - d * sin(((R - r) / r) * 
        (i * THETA));

    }
}

void hypocycloid(double xValues[], double yValues[], int numPoints, 
double R, double r){
    for (int i = 0; i < numPoints; i++) {
        xValues[i] = (R - r) * cos(i * THETA) + r * cos(((R - r) / r) *
        (i * THETA));
        yValues[i] = (R - r) * sin(i * THETA) - r * sin(((R - r) / r) *
        (i * THETA));
    }
}

void butterfly(double xValues[], double yValues[], int numPoints){
    for (int i = 0; i < numPoints; i++) {
        xValues[i] = sin(i * THETA) * ((exp (cos(i * THETA))) - 
        2 * cos(4 * (i * THETA) ) - (pow (sin((i * THETA) / 12), 5)));
        yValues[i] = cos(i * THETA) * ((exp (cos(i * THETA))) - 
        2 * cos(4 * (i * THETA) ) - (pow (sin((i * THETA) / 12), 5)));
    }
}

void writeData(double xValues[], double yValues[], int numPoints,
string pngFilename){
    string filename = "data.dat";

    ofstream dataFile;
    dataFile.open(filename.c_str());

    for (int i = 0; i < numPoints; i++) {
        dataFile << xValues[i] << " " << yValues[i] << endl;
    }

    dataFile.close();

    string plot = "graph.gnuplot";

    ofstream gnuIchooseYou;
    gnuIchooseYou.open(plot.c_str());

    gnuIchooseYou << "set terminal png medium size 640,640 background '#FFFFFF'"    << endl;
    gnuIchooseYou << "set xlabel \"x\"" << endl;
    gnuIchooseYou << "set ylabel \"y\"" << endl;
    gnuIchooseYou << "set key off" << endl;
    gnuIchooseYou << "set output \"tmp.png\"" << endl;
    gnuIchooseYou << "plot \"data.dat\"" << " with lines" << endl;

    gnuIchooseYou.close();

    system("gnuplot graph.gnuplot");
    rename("tmp.png", pngFilename.c_str());
    string write = "display " + pngFilename;
    system(write.c_str());
}

