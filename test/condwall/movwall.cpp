//__INSERT_LICENSE__
// $Id: movwall.cpp,v 1.5 2005/03/31 23:48:37 mstorti Exp $

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
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mov_wall::time_step_pre(double time,int step) { 
  double 
    T = 1,
    Ly = 1,
    Lslit = 0.5*Ly;
  for (int j=0; j<nelem; j++) {
    double y = xwall.e(j,1);
    double R;
    if (y>Lslit) R=1;
    else {
      double y1 = fmod(time,T)/T*Ly;
      double y0 = y1 - Lslit;
      R = (y>y0 && y<y1? 1.0 : 0.0);
    }
    data_p->Rv.ref(j) = R;
    // printf("sliding wall node j %d, y %f, R %f\n",j,y,R);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mov_wall::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mov_wall::close() {}

DL_GENERIC_HOOK(mov_wall);
