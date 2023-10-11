#include <iostream>

using namespace std;

int main ()
{
    // local variable declaration
    int a = 10;

    // do loop execution
    do
    {
        cout << "value of a: " << a << endl;
        a = a + 1;
    }
    while( a < 20 );

    return 0;
}

/* Unlike for and while loops, which test the loop condition at the top of the
 * loop, the do..while loop checks its condition at the bottom of the loop.
 *
 * A do..while loop is similiar to a while loop, except that a do..while loop
 * is guaranteed to execute at least one time.
