#include "Constraint.h"

#include <fem.h>

PyPF::Constraint::Constraint()
  : Ptr(new ::Constraint)
{ 

}

PyPF::Constraint::~Constraint()
{ 
  /* delete[] this->ptr; */
}

