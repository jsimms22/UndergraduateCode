//Name: Jeremiah Simmons
//Description: Program to Calculate Roots of a Number
//Assignment: CSCI-2000 Assignment 3

#include <iostream>
#include <cmath>

using namespace std;

int main () 
{
    cout << "Welcome to the roto-rooter! It's math!" << endl;
    cout << endl;

    double num;  //the number you want to find the roots of
    
    do 
    {
        cout << "Enter a number larger than 0: ";
        cin >> num;
        if (num <= 0) 
        {
            cout << "     Please make it larger than 0" << endl;
        }
    } while (num <= 0);
    
    cout << endl;

    int g;      /*number of iterations the program will 
                run to try to find roots of num*/

    do 
    {
        cout << "Enter maximum tries: ";
        cin >> g;
        if (g <= 0) 
        {
            cout << "     Please make it larger than 0" << endl;
        }
    } while (g <= 0);

    cout << endl;

    int g2 = g - g;             //numbered iteration step for single square root guess
    double newNum2 = num / 2;   //needed to track the guesses
    double oldNum2;             //needed to track the previous guess

    cout << "Square Root" << endl;
    
    while (g2 <= g && abs(newNum2 - oldNum2) >= .000001) 
    {
        oldNum2 = newNum2;
        newNum2 = ((num/oldNum2) + oldNum2) / 2;
        cout << "       Guess " << g2 << ": " << newNum2 << endl;
        g2++;
    }

    cout << "       SQRT(" << num << ") = " << newNum2 << endl;
    cout << endl;

    int g3 = g - g;             //numbered iteration step for single cube root guess
    double newNum3 = num / 2;   //needed to track the guesses
    double oldNum3;             //needed to track the previous guess

    cout << "Cube Root" << endl;

    while (g3 <= g && abs(newNum3 - oldNum3) >= .000001)
    {
        oldNum3 = newNum3;
        newNum3 = ((num / (oldNum3 * oldNum3)) + 2 * oldNum3) / 3;
        cout << "       Guess " << g3 << ": " << newNum3 << endl;
        g3++;
    }

    cout << "       CUBERT(" << num << ") = " << newNum3 << endl;
    cout << endl;

    int g4 = g - g;             //iteration step for single fourth root guess
    double newNum4 = num / 2;   //needed to track the guesses
    double oldNum4;             //needed to track the previous guess

    cout << "Fourth Root" << endl;

    while (g4 <= g && abs(newNum4 - oldNum4) >= .000001)
    {
        oldNum4 = newNum4;
        newNum4 = ((num / (oldNum4 * oldNum4 * oldNum4)) + 3 * oldNum4) / 4;
        cout << "       Guess " << g4 << ": " << newNum4 << endl;
        g4++;
    }

    cout << "       FOURTHRT(" << num << ") = " << newNum4 << endl;
    
    return 0;
}   
