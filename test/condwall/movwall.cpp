//__INSERT_LICENSE__
// $Id: movwall.cpp,v 1.9 2005/04/01 16:08:14 mstorti Exp $

#include <cstdio>
#include <cassert>

#include <map>

#include <src/texthf.h>
#include <src/fem.h>
#include <src/util3.h>
#include <src/hook.h>
#include <src/dlhook.h>
#include <src/dvector.h>
#include "../../applications/ns/condwall.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class mov_wall {
private:
  dvector<double> xwall;
  double Uwall;
  int nelem, ndim;
  cond_wall_data_t *data_p;
public:
  mov_wall() : data_p(NULL) { }
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mov_wall::init(Mesh &mesh_a,Dofmap &dofmap,
	  TextHashTableFilter *options,const char *name) {
  int ierr;
  TGETOPTDEF_ND(GLOBAL_OPTIONS,int,ndim,0);
  assert(ndim>0);

  xwall.cat("condwall.wall-nod.tmp");
  assert(xwall.size() % ndim ==0);
  nelem = xwall.size()/ndim;
  xwall.reshape(2,nelem,ndim);

  string elemset_name = "wall";
  cond_wall_data_map[elemset_name] = cond_wall_data_t();
  data_p = &cond_wall_data_map[elemset_name];
  data_p->Rv.resize(nelem);
  data_p->u1.a_resize(2,nelem,ndim);
  data_p->u2.a_resize(2,nelem,ndim);
  Uwall = 0.2;			// Velocity in `y' direction
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mov_wall::time_step_pre(double time,int step) { 
  double 
    Ly = 1,
    T = Ly/Uwall,
    Lslit = 0.5*Ly;
  for (int j=0; j<nelem; j++) {
    double y = xwall.e(j,1);
    double R = 0;
    // plate1 == fixed plate
    // plate2 == moving plate
    int in_plate1=y>Lslit;

    // Start position of moving plate
    double y0 = Lslit+Uwall*time;
    double dy = modulo(y-y0,Ly);
    int in_plate2 = dy<Lslit;

    // Resistance
    if (in_plate1 || in_plate2) R=1;
    data_p->Rv.ref(j) = R;

    // Velocity on outlet side
    double v=0;
    if (in_plate2) v=Uwall;
    if (in_plate1) v=0.;
    data_p->u2.e(j,1) = v;

    // Velocity on inlet side
    v=0;
    if (in_plate1) v=0.;
    if (in_plate2) v=Uwall;
    data_p->u1.e(j,1) = v;

  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mov_wall::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mov_wall::close() {}

DL_GENERIC_HOOK(mov_wall);
