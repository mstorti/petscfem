// -*- c++ -*-

#ifndef PYPF_FIXATION_H
#define PYPF_FIXATION_H

#include "petscfem4py.h"

PYPF_NAMESPACE_BEGIN

class Fixation
{
 public:
  int* node;
  int* field;
  
  Fixation();
  ~Fixation();

};


PYPF_NAMESPACE_END

#endif // PYPF_FIXATION_H
