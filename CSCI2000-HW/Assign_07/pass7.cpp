//Name: Jeremiah Simmons
//Description: Game of Pig
//Assignment: CSCI - 2000 Assignment 7

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "dice.h"
#include "player.h"

using namespace std;

int main() {
    //Player p1();
    //Player p2();

    //string *player1 = new string(*name);
    //string *player2 = new string((*Player).name);
    string *player1 = new string();
    string *player2 = new string();
    cout << "What are the players' names: \n";
    cin >> *player1 >> *player2;
    Player p1(*player1);
    Player p2(*player2);

    int choice;
    
    cout << "Welcome to Jeremiah's game of Pig\n"
         /*<< "pointers suck\n"*/
         << "Please indicate your choice of sided dye:\n";
    cin >> choice;
    Dice dice(choice);

    int value = 0;
    string i = "no";

    p1.setScore(value);
    p2.setScore(value);
    string answer;

    cout << endl;
    
    for (int t = 0; p1.getScore() < 100 && p2.getScore() < 100; t++) {
        for (int j = 0;i != "yes" && j < 2; j++){
            cout << p1 << "'s roll" << endl;
            cout << "Are you ready to roll?\n";
            cin >> answer;
            while (answer != "yes") {
                cout << "Are you ready to roll now?\n";
                cin >> answer;
            }
            dice.roll();
            value = dice.getValue();
            cout << dice;
            if (j != 1) {
                cout << "Are you happy with your roll?\n";
                cin >> i;
            }
        }
        i = "yes";
        if (i == "yes") {
            p1.setScore(value);
        }
        i = "no";
        cout << "Score is now: " << p1.getScore() << " to " << p2.getScore() 
             << endl;

        for(int j = 0;i != "yes" && j < 2; j++)  {
            cout << p2 << "'s roll" << endl;
            cout << "Are you ready to roll?\n";
            cin >> answer;
            while (answer != "yes") {
                cout << "Are you ready to roll now?\n";
                cin >> answer;
            }
            dice.roll();
            value = dice.getValue();
            cout << dice;
            if (j != 1) {
                cout << "Are you happy with your roll?\n";
                cin >> i;
            }
        }
        i = "yes";
        if(i == "yes") {
            p2.setScore(value);
        }
        i = "no";
        cout << "Score is now: " << p1.getScore() << " to " << p2.getScore()
             << endl;
    }
    
    if (p1.getScore() > p2.getScore()) {
        cout << p1.getName << " is the winner!" << endl;
    } else {
        cout << p2.getName << " is the winner!" << endl;
    }

    delete player1;
    delete player2;
    
}

