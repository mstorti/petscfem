//__INSERT_LICENSE__
//$Id: user.cpp,v 1.1 2002/08/27 16:17:42 mstorti Exp $

#include <iostream>
#include "stack.h"

Stack<char> sc(2);

int main()
{
  Stack<int> si(2);

  try {
    sc.push('a'); sc.push('b');
    si.push(5);   si.push(3);
    cout << "Read back from stack: " << sc.pop() <<
    si.pop() << si.pop() << sc.pop() << endl;

    cout << "Read beyond " << sc.pop() << endl;
  } catch (underflow_error) { cout << "Underflow !\n";
  }


}

