//__INSERT_LICENSE__
//$Id: gasflow.cpp,v 1.3 2003/01/26 21:12:43 mstorti Exp $
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

class gasflow_hook;
typedef map<string,gasflow_hook *> gasflow_hook_table_t;
gasflow_hook_table_t  gasflow_hook_table;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class gasflow_hook {
private:
  double flow_rate,flow_coef;
  string outlet_name;
  int gather_pos;
public:
  double pressure;
  void init(Mesh &mesh,Dofmap &dofmap,
	    TextHashTableFilter *o,const char *name) { 
    TextHashTable *options = (TextHashTable *)TextHashTable::find(name);
    assert(options);
    PetscPrintf(PETSC_COMM_WORLD,
		"-- gasflow_hook::init() name: %s \n",name);
    options->print("options table: ");
    int ierr;
    TGETOPTDEF_ND(options,double,flow_coef,0.);
    TGETOPTDEF(options,double,initial_pressure,0.);
    assert(initial_pressure>0.);
    pressure = initial_pressure;
    TGETOPTDEF_ND(options,int,gather_pos,-1);
    assert(gather_pos>=0);
    TGETOPTDEF_ND(options,double,flow_rate,0.);

    outlet_name = string(name);
    PetscPrintf(PETSC_COMM_WORLD,
		"registering gasflow_hook: %p\n",this);
    assert(gasflow_hook_table.find(outlet_name)
	   ==gasflow_hook_table.end());
    gasflow_hook_table[outlet_name] = this;
  }
  void time_step_pre(double time,int step) { }
  void time_step_post(double time,int step,
		      const vector<double> &gather_values) { 
    double flow_rate_now = gather_values[gather_pos];
    double new_pressure = pressure - flow_coef*(flow_rate_now-flow_rate);
    PetscPrintf(PETSC_COMM_WORLD,
		"flow %g, desired f. %g, p %g, new_p %g\n",
		flow_rate_now,flow_rate,pressure,new_pressure);
    pressure = new_pressure;
  }
  void close() { }
};

DL_GENERIC_HOOK(gasflow_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class flow_controller : public DLGenericTmpl {
private:
  string outlet_name_m;
  gasflow_hook *my_gasflow_hook;
public:
  flow_controller() : my_gasflow_hook(NULL) {}
  void init(TextHashTable *thash) { 
      int ierr;
      TGETOPTDEF_S(thash,string,outlet_name,<none>);
      outlet_name_m = outlet_name;
  }
  double eval(double) { 
    if (!my_gasflow_hook) {
      PetscPrintf(PETSC_COMM_WORLD,
		  " -- flow_controller::init() -- name %s\n",
		  outlet_name_m.c_str());
      gasflow_hook_table_t::iterator q = 
	gasflow_hook_table.find(outlet_name_m);
      assert(q!=gasflow_hook_table.end());
      my_gasflow_hook = q->second;
      assert(my_gasflow_hook);
      PetscPrintf(PETSC_COMM_WORLD,
		  "Looking for hook, gets pointer %p\n",my_gasflow_hook);
    }
    // printf("flow_controller.eval(): returns p=%g\n",
    // my_gasflow_hook->pressure);
    return my_gasflow_hook->pressure;  
  }
};

DEFINE_EXTENDED_AMPLITUDE_FUNCTION2(flow_controller);
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
