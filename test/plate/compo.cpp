//__INSERT_LICENSE__
//$Id: compo.cpp,v 1.1 2002/12/29 21:18:57 mstorti Exp $
#define _GNU_SOURCE

extern int MY_RANK,SIZE;

#include <src/ampli.h>
#include <cstdio>
#include <sys/stat.h>
#include <cassert>
#include <cstdlib>
#include <vector>

#include <src/vecmacros.h>
#include <src/fstack.h>
#include <src/texthash.h>
#include <src/fem.h>
#include "../ROSI/fifo.h"

const int ND=3;

// If needed, this should be set as a class and each communicator
// whould have its own instance.

class composed_movement {
private:
  double Omega,Omega_tr;
  int field;
public:
  composed_movement()  {}
  void init(TextHashTable *thash);
  double eval(double);
  ~composed_movement();
};

void composed_movement::init(TextHashTable *thash) {
  int ierr;
  TGETOPTDEF_ND(thash,double,Omega,0);
  TGETOPTDEF_ND(thash,double,Omega_tr,0.);
  TGETOPTDEF_ND(thash,int,field,0);
  assert(field>=1 && field<=6);
}

double composed_movement::eval(double t) {
  double retval;
  if (field==1)      retval = +sin(Omega_tr*t)*cos(Omega*t);
  else if (field==2) retval = -sin(Omega_tr*t)*sin(Omega*t);
  else if (field==6) retval = Omega;
  else if (field>=1 && field<=6) retval = 0.;
  else {
    PetscPrintf(PETSC_COMM_WORLD,
		"composed_movement: bad field: %d\n",field);
    PetscFinalize();
    exit(0);
  }
//    printf("composed_movement: field %d, t %f, val %f\n",
//  	 field,t,retval);
  return retval;
}

composed_movement::~composed_movement() {}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION(composed_movement);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

class coning {
private:
  double Omega_spin, Omega_coning, coning_angle, gravity;
  int field;
public:
  coning()  {}
  void init(TextHashTable *thash);
  double eval(double);
  ~coning();
};

void coning::init(TextHashTable *thash) {
  int ierr;
  TGETOPTDEF_ND(thash,double,Omega_spin,0);
  TGETOPTDEF_ND(thash,double,Omega_coning,0.);
  TGETOPTDEF_ND(thash,double,coning_angle,0);
  TGETOPTDEF_ND(thash,double,gravity,0);
  TGETOPTDEF_ND(thash,int,field,0);
  assert(field>=1 && field<=6);
}

double coning::eval(double t) {
  double retval;
  double Omega_par,Omega_perp;
  Omega_par = Omega_spin + Omega_coning * cos(coning_angle);
  Omega_perp = Omega_coning * sin(coning_angle);
  if      (field==1) retval = -gravity*sin(coning_angle)*cos(Omega_spin*t);
  else if (field==2) retval = +gravity*sin(coning_angle)*sin(Omega_spin*t);
  else if (field==3) retval = gravity*cos(coning_angle);
  else if (field==4) retval = -Omega_perp*cos(Omega_spin*t);
  else if (field==5) retval = +Omega_perp*sin(Omega_spin*t);
  else if (field==6) retval = Omega_par;
  else 
    if (field>=1 && field<=6) retval = 0.;
  else {
    PetscPrintf(PETSC_COMM_WORLD,
		"coning: bad field: %d\n",field);
    PetscFinalize();
    exit(0);
  }
//    printf("coning: field %d, t %f, val %f\n",
//  	 field,t,retval);
  return retval;
}

coning::~coning() {}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION(coning);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class rot_trans {
private:
  double Omega_spin, Omega_trans, translation_ampl, gravity;
  int field;
public:
  rot_trans()  {}
  void init(TextHashTable *thash);
  double eval(double);
  ~rot_trans();
};

void rot_trans::init(TextHashTable *thash) {
  int ierr;
  TGETOPTDEF_ND(thash,double,Omega_spin,0);
  TGETOPTDEF_ND(thash,double,Omega_trans,0.);
  TGETOPTDEF_ND(thash,double,translation_ampl,0);
  TGETOPTDEF_ND(thash,double,gravity,0);
  TGETOPTDEF_ND(thash,int,field,0);
  assert(field>=1 && field<=6);
}

double rot_trans::eval(double t) {
  double retval;
  if      (field==1) retval = +translation_ampl*cos(Omega_spin*t)*sin(Omega_trans*t);
  else if (field==2) retval = -translation_ampl*sin(Omega_spin*t)*sin(Omega_trans*t);
  else if (field==6) retval = Omega_spin;
  else if (field>=1 && field<=6) retval = 0.;
  else {
    PetscPrintf(PETSC_COMM_WORLD,
		"rot_trans: bad field: %d\n",field);
    PetscFinalize();
    exit(0);
  }
  return retval;
}

rot_trans::~rot_trans() {}

DEFINE_EXTENDED_AMPLITUDE_FUNCTION(rot_trans);
