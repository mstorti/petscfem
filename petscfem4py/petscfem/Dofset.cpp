// $Id: Dofset.cpp,v 1.1.2.6 2006/06/06 16:53:02 dalcinl Exp $

#include "Dofset.h"

#include <fem.h>
#include <dofmap.h>


PYPF_NAMESPACE_BEGIN

Dofset::~Dofset()
{ 
  AmplitudeSet::iterator s = this->amplitude.begin();
  while (s != this->amplitude.end()) { 
    Amplitude* amp = *s++; PYPF_DECREF(amp);
  }
}

Dofset::Dofset()
  : Object(),
    nnod(0), ndof(0),
    amplitude(),
    fixations(),
    constraints()
{ }


Dofset::Dofset(const Dofset& dofset)
  : Object(dofset),
    nnod(dofset.nnod), ndof(dofset.ndof),
    amplitude(dofset.amplitude),
    fixations(dofset.fixations),
    constraints(dofset.constraints)
{ 
  AmplitudeSet::const_iterator s = this->amplitude.begin();
  while (s != this->amplitude.end()) {
    const Amplitude* amp = *s++; PYPF_INCREF(amp);
  }
}

Dofset::Dofset(int nnod, int ndof)
  :  Object(),
     nnod(nnod), ndof(ndof),
     amplitude(),
     fixations(),
     constraints()
{ 
  PYPF_ASSERT(this->nnod>=0,  "invalid 'nnod', out of range (nnod<0)");
  PYPF_ASSERT(this->ndof>=0,  "invalid 'ndof', out of range (ndof<0)");
}

void
Dofset::getSizes(int* nnod, int* ndof) const
{
  if (nnod) *nnod = this->nnod;
  if (ndof) *ndof = this->ndof;
}

void
Dofset::addFixations(int n, 
		     const int    node[], 
		     const int    field[],
		     const double value[])
{
  this->addFixations(n, node, field, value, NULL);
}

void
Dofset::addFixations(int n, 
		     const int    node[], 
		     const int    field[],
		     const double value[], 
		     Amplitude* amplitude)
{
  if (n == 0) return;
  int nnod = this->nnod;
  int ndof = this->ndof;
  for (int i=0; i<n; i++) {
    PYPF_ASSERT(node[i]>=0,    "invalid node, out of range (node<0)");
    PYPF_ASSERT(node[i]<nnod,  "invalid node, out of range (node>=nnod)");
    PYPF_ASSERT(field[i]>=0,   "invalid field, out of range (field<0)");
    PYPF_ASSERT(field[i]<ndof, "invalid field, out of range (field>=ndof)");
  }

  Amplitude::Base* amp = NULL;
  if (amplitude != NULL) {
    if (this->amplitude.insert(amplitude).second) PYPF_INCREF(amplitude);
    amp = *amplitude;
  }
  this->fixations.reserve(this->fixations.size() + n);
  for (int i=0; i<n; i++) {
    this->fixations.push_back(Fixation(node[i], field[i], value[i], amp));
  }
}


void
Dofset::addConstraints(int n, 
		       const int    node[],
		       const int    field[],
		       const double coeff[])
{
  if (n == 0) return;

  PYPF_ASSERT(n>=2, "invalid constraint size");

  int nnod = this->nnod;
  int ndof = this->ndof;
  for (int i=0; i<n; i++) {
    PYPF_ASSERT(node[i]>=0,    "invalid node, out of range (node<0)");
    PYPF_ASSERT(node[i]<nnod,  "invalid node, out of range (node>=nnod)");
    PYPF_ASSERT(field[i]>=0,   "invalid field, out of range (field<0)");
    PYPF_ASSERT(field[i]<ndof, "invalid field, out of range (field>=ndof)");
  }

  this->constraints.push_back(Constraint());
  Constraint& c = this->constraints.back();
  c.reserve(n);
  for (int i=0; i<n; i++) {
    c.push_back(constraint(node[i], field[i], coeff[i]));
  }
}

void
Dofset::view() const {
  printf("Dofset Object:\n");
  printf("  nnod=%d, ndof=%d", this->nnod, this->ndof);
  printf("  Fixations:\n");
  FixationList::const_iterator f = this->fixations.begin();
  while (f != this->fixations.end()) {
    printf("    %d, %d, %g", f->node, f->field, f->value); 
    if (f->amp) printf(", amp: %p", (void*)f->amp);
    printf("\n");
    f++;
  }
  printf("  Constraints:\n");
  ConstraintList::const_iterator cl = this->constraints.begin();
  while (cl != this->constraints.end()) {
    Constraint::const_iterator c = cl->begin();
    printf("    ");
    while (c != cl->end()) {
      printf("%d, %d, %g, ", c->node, c->field, c->coeff); 
      c++;
    }
    printf("\n");
    cl++;
  }
}

void
Dofset::clear()
{
  AmplitudeSet::iterator s = this->amplitude.begin();
  while (s != this->amplitude.end()) { 
    Amplitude* amp = *s++; PYPF_DECREF(amp);
  }
  this->amplitude.clear();
  this->fixations.clear();
  this->constraints.clear();
}

PYPF_NAMESPACE_END
