//__INSERT_LICENSE__
// $Id: mmv-force.cpp,v 1.2 2007/02/11 23:45:34 mstorti Exp $
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
#include <src/h5utils.h>
#include <applications/advdif/mmvforce.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
mmv_force_hook_t::mmv_force_hook_t() : coords_buff(NULL) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*> 
void mmv_force_hook_t::init(Mesh &mesh_a,Dofmap &dofmap_a,
		 TextHashTableFilter *options,
		 const char *name) {
  int ierr;
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,ndim,0);
  assert(ndim>0);

  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,nsaverot,0);
  TGETOPTDEF_S_ND(GLOBAL_OPTIONS,string,xaleppatern,"");
  
  coords_buff = mesh_a.nodedata->nodedata;
  nnod = mesh_a.nodedata->nnod;
  nu = mesh_a.nodedata->nu;
  assert(nu==2*ndim);

  coords0.a_resize(2,nnod,2*ndim).set(coords_buff);
  coords.clone(coords0);
  frame=0;

  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,chimera_save,0);
  // Call the hook for the user
  init_user(mesh_a,dofmap_a,options,name);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void mmv_force_hook_t::time_step_pre(double time,int step) {
  // Shift the coords at t[n] from slot 0 to slot 1.
  // Then call the user function to compute the new coords
  // in slot 0
  for (int j=0; j<nnod; j++) {
    for (int k=0; k<ndim; k++) {
      // Shift old coords to slot '1'
      coords.e(j,ndim+k) = coords.e(j,k);
    }
  }

  // Call the user function to compute the new coords
  // Shift the coords at t[n] from slot 0 to slot 1.
  compute_coords(coords0,coords,time,step);
  
  dvector<double> xale;
  xale.a_resize(2,nnod,ndim);
  // coords is nnod x 2*ndim, we have to copy the first ndim
  // cols to a temporary array xale
  for (int j=0; j<nnod; j++) {
    for (int k=0; k<ndim; k++) {
      xale.e(j,k) = coords.e(j,k);
    }
  }
  if (chimera_save) {
    const char *fname = "./coords.h5";
    if (!access(fname,F_OK)) unlink(fname);
    h5_dvector_write(xale,fname,"xale");
    // FIXME:= the computation of the interpolators should
    // go to the Chimera module
    system("octave-cli -qH mkinterpolators.m > mkinterpolators.log");
  }
  if (xaleppatern.size()>0
      && step%nsaverot==0) {
    // Save coords to H5 file
    // The name of the H5 file for this coords
    char line[1000];
    sprintf(line,xaleppatern.c_str(),frame++);
    h5_dvector_write(xale,line,"xale");
  }
  coords.export_vals(coords_buff);
  // Call the hook for the user
  time_step_pre_user(time,step);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void mmv_force_hook_t::
time_step_post(double time,int step,
	       const vector<double> &gather_values) {
  // Call the hook for the user
  time_step_post_user(time,step,gather_values);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void mmv_force_hook_t::close() {
  // Call the hook for the user
  close_user();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void init_hooks() {
  mmv_force_hook_t MMV_FORCE_HOOK_DUMMY;
}
