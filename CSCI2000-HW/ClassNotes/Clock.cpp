#include "Clock.h"
#include <iostream>

using namespace std;

Clock::Clock() { mValue = 0;}
Clock::Clock(int value) { mValue = value % 12;}
int Clock::getValue() const { return mValue;}
void Clock::setValue(int value) { mValue = value % 12;}

ostream& operator<<(ostream& out, const Clock& cn) //use const as saftey
{   out << cn.getValue() << "cn";                 //want to output not change
    return out;// use out not cout b/c use it with files and not just cout
}
istream& operator>>(istream& in, Clock& cn)
{   int value;
    in >> value;
    cn.setValue(value);
    return in;
}
Clock operator+(Clock &left, int right) 
// more efficient to pass class by reference
{   Clock result(left.getValue() + right);
    return result;
}
Clock operator+(int left, Clock &right) 
{   Clock result(left + right.getValue());
    return result;
}
