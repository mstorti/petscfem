//__INSERT_LICENSE__
// $Id: mmv-force.cpp,v 1.9 2007/02/04 17:00:30 mstorti Exp $
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
#include <src/debug.h>

extern int MY_RANK,SIZE;

#define LINEAR_MMV 1
#define SINE_MMV 2

// #define MMV_FUN LINEAR_MMV
#define MMV_FUN SINE_MMV

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
class mmv_force : public Hook {
private:
  double *coords_buff;
  dvector<double> coords0,coords;
  int ndim, nnod, nu;
  double U,Ly;
public:
  mmv_force() : coords_buff(NULL) { }
  void init(Mesh &mesh_a,Dofmap &dofmap_a,
	    TextHashTableFilter *options,
	    const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close() {}
};


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*> 
void mmv_force::init(Mesh &mesh_a,Dofmap &dofmap_a,
		 TextHashTableFilter *options,
		 const char *name) {
  int ierr;
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,ndim,0);
  assert(ndim>0);

  TGETOPTDEF_ND(GLOBAL_OPTIONS,double,U,NAN);
  assert(!isnan(U));

  TGETOPTDEF_ND(GLOBAL_OPTIONS,double,Ly,NAN);
  assert(!isnan(Ly));

  coords_buff = mesh_a.nodedata->nodedata;
  nnod = mesh_a.nodedata->nnod;
  nu = mesh_a.nodedata->nu;
  assert(nu==3*ndim);

  coords0.a_resize(2,nnod,3*ndim).set(coords_buff);
  coords.clone(coords0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void mmv_force::time_step_pre(double time,int step) {
  double 
    T = 1, 
    omega = 2.0*M_PI/T,
    fac_max = 5.0,
    fac = 1+0.5*(fac_max-1.0)*(1.0-cos(omega*time));
  printf("mmv_force::time_step_pre: time %f, fac %f\n",time,fac);
  for (int j=0; j<nnod; j++) {
    for (int k=0; k<ndim; k++) {
      // Copy old coords in slot `1' to slot `2'
      coords.e(j,2*ndim+k) = coords.e(j,ndim+k);
      // Copy old coords in slot `0' to slot `1'
      coords.e(j,ndim+k) = coords.e(j,k);
    }

    // Set new coords
#if 0
    coords.e(j,0) = 0.5*Ly + fac*(coords0.e(j,0)-0.5*Ly);
    coords.e(j,1) = 0.5*Ly + fac*(coords0.e(j,1)-0.5*Ly);
#else
    double 
      coef = (fac_max-1.0)/(fac_max+1.0),
      xx = coords0.e(j,0),
      yy = coords0.e(j,1),
      tfac = time/(1.0+time),
      shiftx = coef*xx*(1.0-xx)*cos(omega*time)*tfac,
      shifty = coef*yy*(1.0-yy)*cos(omega*time)*tfac;
    coords.e(j,0) = xx + shiftx;
    coords.e(j,1) = yy + shifty;
#endif
  }
  char line[1000];
  sprintf(line,"STEPS/gascont-mmv.state-%d.tmp",step);
  coords.print(line);
  // coords.print();
  coords.export_vals(coords_buff);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void mmv_force::
time_step_post(double time,int step,
	       const vector<double> &gather_values) {
}

DL_GENERIC_HOOK(mmv_force);
