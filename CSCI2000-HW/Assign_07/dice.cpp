#include "dice.h"
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>

using namespace std;

Dice::Dice (int value){ numSides = value; }
int Dice::getSides() const {return numSides;}
int Dice::getValue() const {return topValue;}
void Dice::roll() { srand(time(0)); topValue = rand() % numSides + 1; }

ostream& operator<<(ostream& out, const Dice& dye) {
    if (dye.getSides() != 6) {
    out << "Dice( " << dye.getSides() << ", " << dye.getValue() << ")" << endl;
    } else {
        if (dye.getValue() == 1) {
        out << "+---+\n" << "|   |\n" << "| * |\n" << "|   |\n" << "+---+\n";
        }
        if (dye.getValue() == 2) {
        out << "+---+\n" << "|   |\n" << "|* *|\n" << "|   |\n" << "+---+\n";
        }
        if (dye.getValue() == 3) {
        out << "+---+\n" << "| * |\n" << "|   |\n" << "|* *|\n" << "+---+\n";
        }
        if (dye.getValue() == 4) {
        out << "+---+\n" << "|* *|\n" << "|   |\n" << "|* *|\n" << "+---+\n";
        }
        if (dye.getValue() == 5) {
        out << "+---+\n" << "|* *|\n" << "| * |\n" << "|* *|\n" << "+---+\n";
        }
        if (dye.getValue() == 6) {
        out << "+---+\n" << "|* *|\n" << "|* *|\n" << "|* *|\n" << "+---+\n";
        }
    }
    return out;
}
