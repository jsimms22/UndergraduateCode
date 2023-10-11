#ifndef ClockNumH
#define ClockNumH

#include <iostream>

class Clock
{   private:
        int mValue;
    public:
        Clock();
        Clock(int);
        int getValue() const;
        void setValue(int);

};
Clock operator+(Clock &left, int right);
Clock operator+(int left, Clock &right);
std::ostream& operator<<(std::ostream&, const Clock&);
std::istream& operator>>(std::istream&, Clock&);
#endif

