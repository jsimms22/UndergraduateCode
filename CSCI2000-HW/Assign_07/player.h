#ifndef Person_H
#define Person_H

#include <iostream>
#include <string>

using namespace std;

class Player {
    private:
        string mName;
        int score;
    public:
        Player (string);
        string getName() const;
        int getScore() const;
        void setScore(int);
};
ostream& operator<<(ostream&, const Player&);
#endif
