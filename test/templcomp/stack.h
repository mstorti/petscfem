#ifndef STACK_H
#define STACK_H

#include <stdexcept>
// create a stack class for type T
#pragma interface
// a generic stack class

template<class T> class Stack {
  T* v;          // stack array
  int max_size;  // size of stack array
  int top;       // points to first free element
 public:
  class underflow_error {  };
  class overflow_error {  };

  Stack (int size); // constructor for Stack of size 'size'
  ~Stack();         // destructor

  void push(T);
  T pop();
};

typedef Stack<int> Stackint;
typedef Stack<char> Stackchr;


#endif
