#include "player.h"
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

Player::Player (string name) { mName = name; score = 0; }
string Player::getName() const { return mName; }
int Player::getScore() const { return score; }
void Player::setScore(int value) { score = score + value; }

ostream& operator<<(ostream& out, const Player& p1) {
    out << p1.getName();
    return out;
}


