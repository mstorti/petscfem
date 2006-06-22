// $Id: Dofset.cpp,v 1.1.2.9 2006/06/22 22:34:52 dalcinl Exp $

#include "Dofset.h"

#include <fem.h>
#include <dofmap.h>


PYPF_NAMESPACE_BEGIN

Dofset::~Dofset()
{ }

Dofset::Dofset()
  : Object(),
    nnod(0), ndof(0),
    fixations(),
    amplitudes(),
    constraints()
{ }


Dofset::Dofset(const Dofset& dofset)
  : Object(dofset),
    nnod(dofset.nnod), ndof(dofset.ndof),
    fixations(dofset.fixations),
    amplitudes(dofset.amplitudes),
    constraints(dofset.constraints)
{ }

static void chk_sizes(int nnod, int ndof) {
  PYPF_ASSERT(nnod>=0, "invalid 'nnod', out of range (nnod<0)");
  PYPF_ASSERT(ndof>=0, "invalid 'ndof', out of range (ndof<0)");
}

Dofset::Dofset(int nnod, int ndof)
  :  Object(),
     nnod(nnod), ndof(ndof),
     fixations(),
     amplitudes(),
     constraints()
{ chk_sizes(this->nnod, this->ndof); }

Dofset::Dofset(int nnod, int ndof, MPI_Comm comm)
  :  Object(comm),
     nnod(nnod), ndof(ndof),
     fixations(),
     amplitudes(),
     constraints()
{ chk_sizes(this->nnod, this->ndof); }


void
Dofset::chk_fixa(int n, int f)
{
  PYPF_ASSERT(n >= 0,          "invalid node, out of range (node<0)");
  PYPF_ASSERT(n <  this->nnod, "invalid node, out of range (node>=nnod)");
  PYPF_ASSERT(f >= 0,          "invalid field, out of range (field<0)");
  PYPF_ASSERT(f <  this->ndof, "invalid field, out of range (field>=ndof)");
}

void
Dofset::chk_fixa(int n, int f, double v)
{
  this->chk_fixa(n, f);
  PYPF_ASSERT(!(v!=v), "invalid value, not a number (NaN)");
}

void
Dofset::chk_fixa(int n, 
		 const int node[],
		 const int field[]) 
{
  for (int i=0; i<n; i++)
    this->chk_fixa(node[i], field[i]);
}

void
Dofset::chk_fixa(int n, 
		 const int    node[],
		 const int    field[],
		 const double value[]) 
{
  for (int i=0; i<n; i++)  
    this->chk_fixa(node[i], field[i], value[i]);
}

void
Dofset::add_fixa(int n, int f, double v, Amplitude::Base* a)
{
  this->fixations.push_back(Fixation(n, f, v, a));
}

void
Dofset::add_fixa(int n, 
		 const int    node[],
		 const int    field[],
		 const double value[],
		 Amplitude::Base* a)
{
  if (n == 0) return;
  this->fixations.reserve(this->fixations.size() + n);
  for (int i=0; i<n; i++)
    this->add_fixa(node[i], field[i], value[i], a);
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
		     const double value[],
		     Amplitude* amplitude)
{
  Amplitude::Base* amp = NULL;
  if (amplitude) amp = *amplitude;
  this->chk_fixa(n, node, field, value);
  this->add_fixa(n, node, field, value, amp);
  this->amplitudes.add(amplitude);
}

void
Dofset::addConstraints(int n, 
		       const double coeff[],
		       const int    node[],
		       const int    field[])
{
  if (n == 0) return;
  PYPF_ASSERT(n>=2, "invalid constraint size");
  for (int i=0; i<n; i++) {
    PYPF_ASSERT(!(coeff[i]!=coeff[i]),  "invalid coefficient, not a number (NaN)");
    PYPF_ASSERT(node[i]  >= 0,          "invalid node, out of range (node<0)");
    PYPF_ASSERT(node[i]  < this->nnod,  "invalid node, out of range (node>=nnod)");
    PYPF_ASSERT(field[i] >= 0,          "invalid field, out of range (field<0)");
    PYPF_ASSERT(field[i] < this->ndof,  "invalid field, out of range (field>=ndof)");
  }
  this->constraints.push_back(Constraint());
  Constraint& c = this->constraints.back();
  c.reserve(n);
  for (int i=0; i<n; i++) {
    c.push_back(constraint(coeff[i], node[i], field[i]));
  }
}

void
Dofset::view() const {
  printf("Dofset Object:\n");
  printf("  nnod=%d, ndof=%d\n", this->nnod, this->ndof);
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
      printf("%g, %d, %d, ", c->coeff, c->node, c->field);
      c++;
    }
    printf("\n");
    cl++;
  }
}

void
Dofset::clear()
{
  this->fixations.clear();
  this->amplitudes.clear();
  this->constraints.clear();
}

PYPF_NAMESPACE_END
