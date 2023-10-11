//Name: Jeremiah Simmons
//Description: Program to generate psuedo random geometric art
//Assignment: CSCI 2000 Assignment 4

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

//function for color of every shape
void fill(ofstream &output, int red, int green, int blue);

void askBackground(string &backgroundColor);

//function for up and down triangles
void tri(ofstream &output, double size, int i, int o, int l, int shapeCount,
int shapeType, double x1, double x2, double x3, double y1, double y2,
double y3);

//function for rectangles
void rect(ofstream &output, double size, int i, int o, int shapeCount, 
double x, double y, double width, double height);

void circle(ofstream &output, double size, int i, int o, int shapeCount,
double cx, double cy, double radius);

int main () {

    srand(time(0));

    ofstream output;
    output.open("pass4.svg");

    double size;
    cout << "Enter the image size: ";
    cin >> size;

    output << "<svg xmlns='http://www.w3.org/2000/svg' width='" << size << 
    "' height='" << size << "' inversion='1.1'>" << endl;

    string backgroundColor;
    askBackground(backgroundColor);

    output << "<rect x='0' y='0' width='" << size << "' height='" << size << 
    "' fill='" << backgroundColor << "' />" << endl;
    
    int shapeCount; //used for tracking how many rows/columns are desired

    cout << "Enter the number of shapes per row/column: ";
    cin >> shapeCount;

    int shapeType;
    cout << "Enter the shape(s) to be used\n"
    "   1-square\n"
    "   2-circle\n"
    "   3-up triangle\n"
    "   4-down triangle\n"
    "   5-all\n"
    "   6-all w/blanks\n"
    "choice: ";
    cin >> shapeType;

    int red, green, blue; //variables for fill function
    int l; //used for random number variable in tri() function

    //below are various variables for the other shape functions
    double cx, cy, radius, x1, x2, x3, y1, y2, y3, x, y, width, height;
    
    for(int o = 0; o < shapeCount; o++) {
        for(int i = 0; i < shapeCount; i++) {
            if(shapeType <= 4) {
                switch(shapeType) {
                    case 1:
                        rect(output, size, i, o, shapeCount, x, y, width, height
                        ); 
                        fill(output, red, green, blue); 
                    break;
                    case 2:
                        circle(output, size, i, o, shapeCount, cx, cy, radius);
                        fill(output, red, green, blue);
                    break;
                    case 3:
                        tri(output, size, i, o, l, shapeCount, shapeType, x1, x2
                        , x3, y1, y2, y3);
                        fill(output, red, green, blue);
                    break;
                    case 4:
                        tri(output, size, i, o, l, shapeCount, shapeType, x1, x2
                        , x3, y1, y2, y3);
                        fill(output, red, green, blue);
                    break;
                }
            } if(shapeType == 5) {
                int j = rand() % 3;
                switch(j) {
                    case 0:
                        rect(output, size, i, o, shapeCount, x, y, width, height);
                        fill(output, red, green, blue);
                    break;
                    case 1:
                        circle(output, size, i, o, shapeCount, cx, cy, radius);
                        fill(output, red, green, blue);
                    break;
                    case 2:
                        tri(output, size, i, o, l, shapeCount, shapeType, x1, x2
                        , x3, y1, y2, y3);
                        fill(output, red, green, blue);
                    break;
                    default:
                        tri(output, size, i, o, l, shapeCount, shapeType, x1, x2
                        , x3, y1, y2, y3);
                        fill(output, red, green, blue);
                    break;
                }
            } if(shapeType == 6) {
                int k = rand() % 5;
                switch(k) {
                    case 0:
                        rect(output, size, i, o, shapeCount, x, y, width,
                        height);
                        fill(output, red, green, blue);
                    break;
                    case 1:
                        circle(output, size, i, o, shapeCount, cx, cy, radius
                        );
                        fill(output, red, green, blue);
                    break;
                    case 2:
                        tri(output, size, i, o, l, shapeCount, shapeType, x1, x2
                        , x3, y1, y2, y3);
                        fill(output, red, green, blue);
                    break;
                    case 3:
                        tri(output, size, i, o, l, shapeCount, shapeType, x1, x2
                        , x3, y1, y2, y3);
                        fill(output, red, green, blue);
                    break;
                    default:
                        output << endl;
                    break;
                }
            }
        } 
    }

    output << "</svg>" << endl;

    output.close();
}

void askBackground(string &backgroundColor) {

    cout << "Enter the background color: ";
    cin >> backgroundColor;
    
    cout << "   " << backgroundColor << " = ";

    if(backgroundColor == "white") {
        backgroundColor = "rgb(255,255,255)";
    } if(backgroundColor == "blue") {
        backgroundColor = "rgb(0,0,255)";
    } if(backgroundColor == "red") {
        backgroundColor = "rgb(255,0,0)";
    } if(backgroundColor == "green") {
        backgroundColor = "rgb(0,255,0)";
    } else { 
        backgroundColor = "rgb(0,0,0)";
    }

    cout << backgroundColor << endl;
}

void tri(ofstream &output, double size, int i, int o, int l, int shapeCount,
int shapeType, double x1, double x2, double x3, double y1, double y2, 
double y3) {
    
    if(shapeType == 3) {
        x1 = (size / shapeCount) - ((size / shapeCount) / 2) + (i * (size / 
        shapeCount));
        y1 = (size / shapeCount) * (o);

        x2 = (size / shapeCount) * (i + 1);
        y2 = (size / shapeCount) * (o + 1);

        x3 = (size / shapeCount) * (i);
        y3 = (size / shapeCount) * (o + 1);

        output << "<polygon points='" << x1 << "," << y1 << " " << x2 << ","
        << y2 << " " << x3 << "," << y3 << "' ";
    } if(shapeType == 4) {
        x1 = (size / shapeCount) - ((size / shapeCount) / 2) + (i * (size / 
        shapeCount));
        y1 = (size / shapeCount) * (o + 1);

        x2 = (size / shapeCount) * (i + 1);
        y2 = (size / shapeCount) * (o);

        x3 = (size / shapeCount) * (i);
        y3 = (size / shapeCount) * (o);

        output << "<polygon points='" << x1 << "," << y1 << " " << x2 << ","
        << y2 << " " << x3 << "," << y3 << "' ";
    } if(shapeType == 5 || shapeType == 6) {
        l = rand() % 2;
        if (l == 0){
            x1 = (size / shapeCount) - ((size / shapeCount) / 2) + (i * (size / 
            shapeCount));
            y1 = (size / shapeCount) * (o);

            x2 = (size / shapeCount) * (i + 1);
            y2 = (size / shapeCount) * (o + 1);

            x3 = (size / shapeCount) * (i);
            y3 = (size / shapeCount) * (o + 1);

            output << "<polygon points='" << x1 << "," << y1 << " " << x2 << ","
            << y2 << " " << x3 << "," << y3 << "' ";
        } if (l == 1) {
            x1 = (size / shapeCount) - ((size / shapeCount) / 2) + (i * (size / 
            shapeCount));
            y1 = (size / shapeCount) * (o + 1);

            x2 = (size / shapeCount) * (i + 1);
            y2 = (size / shapeCount) * (o);

            x3 = (size / shapeCount) * (i);
            y3 = (size / shapeCount) * (o);

            output << "<polygon points='" << x1 << "," << y1 << " " << x2 << ","
            << y2 << " " << x3 << "," << y3 << "' ";
        }
    }
}

void fill(ofstream &output, int red, int green, int blue) {
    red = rand() % 256;
    green = rand() % 256;
    blue = rand() % 256;
    output << "fill='rgb(" << red << "," << green << "," << blue <<")' />\n";
}

void circle(ofstream &output, double size, int i, int o, int shapeCount,
double cx, double cy, double radius) {
    cx = (size / shapeCount) - ((size / shapeCount) / 2) + (i * (size / shapeCount));
    cy = (size / shapeCount) - ((size / shapeCount) / 2) + (o * (size / shapeCount));
    radius = (size / shapeCount) / 2;

    output << "<circle cx='" << cx << "' cy='" << cy << "' r='"
    << radius << "' ";
}

void rect(ofstream &output, double size, int i, int o, int shapeCount, 
double x, double y, double width, double height) {
    x = (size / shapeCount) * (i);
    y = (size / shapeCount) * (o);
    width = (size / shapeCount);
    height = (size / shapeCount);

    output << "<rect x='" << x << "' y='" << y << "' width='" << width <<
    "' height='" << height << "' ";
}
