//__INSERT_LICENSE__
// $Id: mesh-move.cpp,v 1.10 2006/05/13 21:14:29 mstorti Exp $
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <cstdio>
#include <cassert>
#include <cmath>

#include <map>

#include <src/texthf.h>
#include <src/fem.h>
#include <src/util3.h>

#include <src/getprop.h>
#include <src/ampli.h>
#include <src/hook.h>
#include <src/dlhook.h>
#include <src/dvecpar.h>
#include <src/penalize.h>
#include <src/debug.h>

extern int MY_RANK,SIZE;

extern Mesh *GLOBAL_MESH;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class fixa : public DLGenericTmpl {
private:
  double A, omega, trelax, Tend, Tstop;
  int fun_type;
public:
  fixa() { }
  void init(TextHashTable *thash);
  double eval(double);
  ~fixa() { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void fixa::
init(TextHashTable *thash) { 
  int ierr;
  //o The amplitud of the shaking at the base
  TGETOPTDEF_ND(thash,double,A,NAN);
  PETSCFEM_ASSERT(A!=NAN && A>=0.0,"bad Amplitude, A= %g",A);

  //o The frequency of the shaking at the base
  TGETOPTDEF_ND(thash,double,omega,NAN);
  PETSCFEM_ASSERT(omega!=NAN && omega>=0.0,
                  "bad frequency, omega= %g",omega);

  //o A relax time to reach the nominal amplitud
  TGETOPTDEF_ND(thash,double,trelax,NAN);
  PETSCFEM_ASSERT(trelax!=NAN && trelax>=0.0,
                  "bad trelax, trelax= %g",trelax);

  TGETOPTDEF_ND(thash,double,Tend,NAN);
  PETSCFEM_ASSERT0(!isnan(Tend),"Tend is required");

  TGETOPTDEF_ND(thash,double,Tstop,NAN);
  PETSCFEM_ASSERT0(!isnan(Tstop),"Tstop is required");

  //o Index for defining type of fun
  TGETOPTDEF_ND(thash,int,fun_type,0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double fixa::eval(double time) { 
  int f = field();
  double dx;
  if (f!=1) dx = 0.0;
  else {
    double w0=1.0, w1=3.0, A0=0.15, A1=A0;
    double 
      tfac = (time<Tend ? time/Tend : 1.0),
      w = w0+tfac*(w1-w0),
      A = A0+tfac*(A1-A0);
    dx = A*sin(w*time);
  }
  return dx;
}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION2(fixa);
