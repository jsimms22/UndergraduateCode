//Name: Jeremiah Simmons
//Description: Program to find roota for polynomials up to 3rd degree
//Assignment: CSI-2000 Assignment 6

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

void evaluateFx(double *a, double *b, double *c, double *d, double *x, double *fxPtr,double *fxDerivPtr) {
    *fxPtr = *a * pow(*x,3) + *b * pow(*x,2) + (*c * *x) + *d;

    *fxDerivPtr = 3 * *a * pow(*x,2) + 2 * (*b * *x) + *c;
}

void nextVal(double *fxPtr, double *fxDerivPtr, double *xPtr, double *xprevPtr){
    *xprevPtr = *xPtr;
    *xPtr = *xPtr - (*fxPtr / *fxDerivPtr);
}

int main () {

    double a, b, c, d;
    
    cout << "Welcome to 'Jeremiah's Root Finder Thingy'\n"
            "Please give us the coefficients of your 3rd order polynominal:\n";
    cout << "a = ";
    cin >> a;
    cout << "b = ";
    cin >> b;
    cout << "c = ";
    cin >> c;
    cout << "d = ";
    cin >> d;

    double *aPtr;
    aPtr = new double(a);
    double *bPtr;
    bPtr = new double(b);
    double *cPtr;
    cPtr = new double(c);
    double *dPtr;
    dPtr = new double(d);

    cout << "Your polynominal is: "
    << a << "x^3 + " << b << "x^2 + " << c << "x + " << d << endl;

    double x;

    cout << "For Newton's method this program requires an initial value for x\n"
            "x = ";
    cin >> x;

    double *xPtr = &x;
    double xprev = 0;
    double *xprevPtr = &xprev;
    double fx, fxDeriv; 
    double *fxPtr;
    fxPtr = &fx;
    double *fxDerivPtr;
    fxDerivPtr = &fxDeriv;
    double diff = 1;
    
    cout << setprecision(5) << "N " 
             << setw(10) << "xN " 
             << setw(10) << "f(xN) " 
             << setw(10) << "f'(xN) "
             << setw(10) << "xN+1 " 
             << setw(10) << "abs(xN - xN+1)" << endl;

    for (int N = 0; diff > .00001 && N<=30;N++) {
        evaluateFx(aPtr,bPtr,cPtr,dPtr,xPtr,fxPtr,fxDerivPtr);
        nextVal(fxPtr,fxDerivPtr,xPtr, xprevPtr);
        diff = abs(*xprevPtr - *xPtr);
        cout << setprecision(5) << N << " " 
             << setw(10) << setprecision(5) << *xprevPtr << " " 
             << setw(10) << setprecision(5) << *fxPtr << " " 
             << setw(10) << setprecision(5) << *fxDerivPtr << " " 
             << setw(10) << setprecision(5) << *xPtr << " " 
             << setw(10) << setprecision(5) << diff << endl;
    }

    cout << "One of your roots is: " << x << endl;

    delete aPtr;
    delete bPtr;
    delete cPtr;
    delete dPtr;
}
