//__INSERT_LICENSE__
//$Id: inviscid.cpp,v 1.1 2002/12/29 21:18:57 mstorti Exp $
#define _GNU_SOURCE

extern int MY_RANK,SIZE;

#include "./inviscid.h"

extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_visc_hook::init(Mesh &mesh,Dofmap &dofmap,
		     TextHashTableFilter *options,const char *name) {
  // File to send forces and moments to "PFM"
  int ierr;

  if (!MY_RANK) {
    printf("VISCOUS: Opening fifos for communicating with INVISCID.\n");

    visc2inv = fopen("visc2inv","w");
    assert(visc2inv);
    setvbuf(visc2inv,NULL,_IOLBF,0);

//      inv2visc = fopen("inv2visc","r");
//      assert(inv2visc);

    printf("PETSCFEM: Done.\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_visc_hook::time_step_pre(double t,int step) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void coupling_visc_hook::time_step_post(double time,int step,
			       const vector<double> &gather_values) {
  fprintf(visc2inv,"Computed step %d, time %f\n",step,time);
}

DL_GENERIC_HOOK(coupling_visc_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
