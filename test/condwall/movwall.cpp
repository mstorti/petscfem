//__INSERT_LICENSE__
// $Id: movwall.cpp,v 1.1 2005/03/29 04:01:36 mstorti Exp $

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
public:
  void init(Mesh &mesh_a,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

void mov_wall::init(Mesh &mesh_a,Dofmap &dofmap,
	  TextHashTableFilter *options,const char *name) {
}

void mov_wall::time_step_pre(double time,int step) { 
  int nelem = cond_wall_resistance.size();
  for (int j=0; j<nelem; j++)
    cond_wall_resistance.ref(j) = (j<nelem/2);
}

void mov_wall::time_step_post(double time,int step,
			      const vector<double> &gather_values) {
}

void mov_wall::close() {}

DL_GENERIC_HOOK(mov_wall);
