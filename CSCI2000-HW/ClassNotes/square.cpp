#include <iostream>
#include <string>

using namespace std;

class Square {
    private:
        int *mSize;
        string *mColor;
    public:
        Square();
        Square(int, string); //Constructor

        ~Square(); //Destructor

        void setSize(int);
        int getSize() const;

        void setColor(string);
        string getColor() const;
};

Square::Square() {
    mSize = new int(1);
    mColor = new string("red");
}

Square::Square(int size, string color) {
    mSize = new int(size);
    mColor = new string(color);
}

Square::~Square() {
    cout << "Destructor called\n";
    delete mSize;
    delete mColor;
}

void Square::setSize(int size) { 
    *mSize = size; 
}

int Square::getSize() const {
    return *mSize;
}

void Square::setColor(string color) { 
    *mColor = color;
}

string Square::getColor() const { 
    return *mColor; 
}

void test() {
    Square sq = Square(10, "blue");

    cout << "The area is " << sq.getSize() * sq.getSize() << endl;
    cout << "The perimeter is " << 4 * sq.getSize() << endl;

    updateSize(sq);
}

void updateSize(Square s) {
    cout << "Enter a new size: ";
    int newSize;
    cin >> newSize;
    s.setSize(newSize);
}

int main() {
    cout << "Better test()\n";
    test();
    cout << "After test()\n";

}
