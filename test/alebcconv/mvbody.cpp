//__INSERT_LICENSE__
// $Id: mvbody.cpp,v 1.2 2006/06/19 12:35:39 mstorti Exp $
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
#include <src/dvector.h>
#include <applications/advdif/advective.h>
#include <src/penalize.h>
#include <src/debug.h>

extern int MY_RANK,SIZE;

// FIXME:= Falta imponer la velocidad sobre el
// contorno del cuerpo

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class mvbodyh : public Hook {
private:
  Mesh *meshp;
  double Dt, vmesh;
  int ndim, nu, nnod;
  Debug debug;
  dvector<double> coords;
public:
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,
	    const char *name);
  void time_step_pre(double time,int step);
  mvbodyh() : debug(0,PETSC_COMM_WORLD) { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mvbodyh
::init(Mesh &mesh_a,Dofmap &dofmap,
       TextHashTableFilter *options,
       const char *name) {
  meshp = &mesh_a;

  int ierr;
#if 0  
  debug.activate();
  Debug::init();
  debug.trace("activate debugging");
#endif

  //o Time step. 
  TGETOPTDEF_ND(GLOBAL_OPTIONS,double,Dt,0);
  assert(Dt>0.0);

  //o Mesh velocity
  TGETOPTDEF_ND(GLOBAL_OPTIONS,double,vmesh,NAN);
  assert(!isnan(vmesh));

  //o Dimension of problem. 
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,ndim,0);
  assert(ndim>0);

  double *xp = meshp->nodedata->nodedata;
  nu = meshp->nodedata->nu;
  nnod = meshp->nodedata->nnod;
  assert(nu==2*ndim);
  coords.mono(nu*nnod)
    .reshape(2,nnod,nu).set(xp);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mvbodyh
::time_step_pre(double time,int step) { 
  double *xp = meshp->nodedata->nodedata;

  for (int j=0; j<nnod; j++) {
    // Copy old coords to ndim+[0,ndim) range
    xp[j*nu+ndim] = xp[j*nu];
    // Compute new coords 
    xp[j*nu] = xp[j*nu+ndim] + vmesh*Dt;
  }
}

DL_GENERIC_HOOK(mvbodyh);
