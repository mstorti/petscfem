//__INSERT_LICENSE__
//$Id: stack.cpp,v 1.1 2002/08/27 16:17:42 mstorti Exp $

// create a stack class for type T
#pragma implementation
#include "stack.h"

template<class T> Stack<T>::Stack(int size)
{
  top = 0;
  max_size = size;
  v = new T[size];
}

// destructor

template<class T> Stack<T>::~Stack()
{
  delete [] v;
}

// push one element on top of stack

template<class T> void Stack<T>::push(T c)
{

  if (top == max_size) throw overflow_error();
  v[top++] = c;
}

// pop one element from top of stack

template<class T> T Stack<T>::pop()
{
  if (top == 0) throw underflow_error();
  return v[--top];
}
