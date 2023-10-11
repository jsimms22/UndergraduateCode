#ifndef Dice_H
#define Dice_H

#include <iostream>

using namespace std;

class Dice {
    private:
        int numSides;
        int topValue;

    public:
        Dice (int);
        int getSides() const;
        int getValue() const;
        void roll();
};
ostream& operator<<(ostream&, const Dice&);
#endif
