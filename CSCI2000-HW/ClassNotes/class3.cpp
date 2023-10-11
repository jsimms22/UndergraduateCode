
//This program is blah

#include <iostream>
#include "Clock.h"
#include <fstream>

using namespace std;

int main()
{
    Clock num(10);
    cout << num.getValue() << endl;
    //num.setValue(num.getValue() + 5);
    num = num + 5;
    cout << num.getValue() << endl;
    num = 7 + num;
    cout << num << endl; // we can do this after operator>>

    cout << "Enter a  clock number: ";
    cin >> num;
    cout << num << endl;
    ofstream dataFiles;
    dataFile.open("data.txt");
    
        for(int i=0; i <6; i++)
        {   datafile << num << endl;
            num = num + 1;
        }

    dataFile.close();
}
