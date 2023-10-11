/***
    Jeremiah Simmons
    CSCI 2000
    Programming assignment 2
***/

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

int main ()
{
    //Program to calculate distance between two coordinates in a 3D space

    double x_initial;   //first coordinate x value
    double y_initial;   //first coordinate y value
    double z_initial;   //first coordinate z value

    double x_final;   //second coordinate x value
    double y_final;   //second coordinate y value
    double z_final;   //second coordinate z value

    cout << "Enter point1 x: ";
    cin >> x_initial;
    cout << "Enter point1 y: ";
    cin >> y_initial;
    cout << "Enter point1 z: ";
    cin >> z_initial;

    cout << endl;

    cout << "Enter point2 x: ";
    cin >> x_final;
    cout << "Enter point2 y: ";
    cin >> y_final;
    cout << "Enter point2 z: ";
    cin >> z_final;

    cout << endl;

    double distance;   

    //Distance formula for two points in a 3D space

    distance = sqrt((x_final - x_initial)*(x_final - x_initial) +
    (y_final - y_initial)*(y_final - y_initial) +
    (z_final - z_initial)*(z_final - z_initial));

    cout << "The distance between ";
    cout << "(" << x_initial << ", " << y_initial << ", " << z_initial << ")";
    cout << " and ";
    cout << "(" << x_final << ", " << y_final << ", " << z_final << ")";
    cout << " is " << setprecision(8) << distance << endl;

    return 0;
}
