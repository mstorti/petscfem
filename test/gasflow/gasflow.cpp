//__INSERT_LICENSE__
//$Id: gasflow.cpp,v 1.8 2003/01/30 18:02:06 mstorti Exp $
#define _GNU_SOURCE

extern int MY_RANK,SIZE;
#include <src/ampli.h>
#include <cstdio>
#include <sys/stat.h>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <map>

#include <src/vecmacros.h>
#include <src/fstack.h>
#include <src/texthash.h>
#include <src/fem.h>
#include <src/texthf.h>

#include <src/hook.h>
#include <src/dlhook.h>


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class cut_regulator_hook;
class flow_controller;

typedef map<string,cut_regulator_hook *> cut_regulator_hook_table_t;
cut_regulator_hook_table_t  cut_regulator_hook_table;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class cut_regulator_hook {
private:
  int nstream;
  double *flow_rate,flow_coef, *flow_rate_now;
  int *gather_pos;
public:
  friend class flow_controller;
  double *pressure;
  cut_regulator_hook() :
    flow_rate(NULL), gather_pos(NULL), pressure(NULL),
    flow_rate_now(NULL) {}
  ~cut_regulator_hook() {
    delete[] flow_rate;
    delete[] flow_rate_now;
    delete[] pressure;
    delete[] gather_pos;
  }
  void init(Mesh &mesh,Dofmap &dofmap,
	    TextHashTableFilter *o,const char *name) { 
    TextHashTable *options = (TextHashTable *)TextHashTable::find(name);
    assert(options);
    PetscPrintf(PETSC_COMM_WORLD,
		"-- cut_regulator_hook::init() name: %s \n",name);
    options->print("options table: ");
    int ierr;
    TGETOPTDEF_ND(options,int,nstream,0);
    assert(nstream>0);

    flow_rate = new double[nstream];
    flow_rate_now = new double[nstream];
    pressure = new double[nstream];
    gather_pos = new int[nstream];

    TGETOPTDEF_ND(options,double,flow_coef,0.);
    ierr = get_double(options,"initial_pressure",pressure,1,nstream); 
    assert(!ierr);
    ierr = get_int(options,"gather_pos",gather_pos,1,nstream); 
    assert(!ierr);
    ierr = get_double(options,"flow_rate",flow_rate,1,nstream); 
    assert(!ierr);
    PetscPrintf(PETSC_COMM_WORLD,
		"registering cut_regulator_hook: %p\n",this);
    string my_name(name);
    assert(cut_regulator_hook_table.find(my_name)
	   ==cut_regulator_hook_table.end());
    cut_regulator_hook_table[my_name] = this;
  }
  void time_step_pre(double time,int step) { }
  void time_step_post(double time,int step,
		      const vector<double> &gather_values) { 
    // double sum_flow_rate_now=0., sum_flow_rate=0., p_avrg=0.;
    for (int j=0; j<nstream; j++) {
      flow_rate_now[j] = (j==0? 1. : -1.) * gather_values[gather_pos[j]];
    }
    assert(flow_rate_now[0]>0.);
    for (int j=1; j<nstream; j++) {
      assert(flow_rate_now[j]<0.);
      double new_pressure = pressure[j] 
	+ flow_coef*(flow_rate_now[j]-flow_rate[j]*flow_rate_now[0]);
      PetscPrintf(PETSC_COMM_WORLD,
		  "stream %d, flow%d = %g, desired flow = %g, p%d = %g, new_p%d = %g\n",
		  j,flow_rate_now[j],j,flow_rate[j]*flow_rate_now[0],
		  j,pressure[j],j,new_pressure);
      pressure[j] = new_pressure;
    }
  }
  void close() { }
};

DL_GENERIC_HOOK(cut_regulator_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class flow_controller : public DLGenericTmpl {
private:
  string cut_regulator_m;
  cut_regulator_hook *my_cut_regulator_hook;
  int index;
public:
  flow_controller() : my_cut_regulator_hook(NULL) {}
  void init(TextHashTable *thash) { 
      int ierr;
      TGETOPTDEF_S(thash,string,cut_regulator,<none>);
      cut_regulator_m = cut_regulator;
      TGETOPTDEF_ND(thash,int,index,-1);
      assert(index>-1);
  }
  double eval(double) { 
    if (!my_cut_regulator_hook) {
      PetscPrintf(PETSC_COMM_WORLD,
		  " -- flow_controller::init() -- name %s\n",
		  cut_regulator_m.c_str());
      cut_regulator_hook_table_t::iterator q = 
	cut_regulator_hook_table.find(cut_regulator_m);
      assert(q!=cut_regulator_hook_table.end());
      my_cut_regulator_hook = q->second;
      assert(my_cut_regulator_hook);
      PetscPrintf(PETSC_COMM_WORLD,
		  "Looking for hook, gets pointer %p\n",my_cut_regulator_hook);
      assert(index < my_cut_regulator_hook->nstream);
    }
    return my_cut_regulator_hook->pressure[index];
  }
};

DEFINE_EXTENDED_AMPLITUDE_FUNCTION2(flow_controller);
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
