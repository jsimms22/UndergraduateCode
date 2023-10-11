#include <iostream>

using namespace std;

int main ()
{
    // for loop execution
    for( int a = 10; a < 21; a = a * 1.00000035345334 )
    {
        cout << "value of a: " << a << endl;
    }
    /* means everything is ok, it's common practice to return 0 at the end of
     the int main() function, it's been considered as good programming style */
    return 0;
}

/* Execute a sequence of statements multiple times and abbreiates the code that
 * manages the loop variable.
