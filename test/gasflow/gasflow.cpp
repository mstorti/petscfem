//__INSERT_LICENSE__
//$Id: gasflow.cpp,v 1.7 2003/01/30 17:01:15 mstorti Exp $
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

#if 0
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
    if (!MY_RANK) options->print("options table: ");
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
    double new_pressure = pressure + flow_coef*(flow_rate_now-flow_rate);
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
#else

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
    // Normalize flow_rate
    double flow_rate_sum=0.;
    for (int j=0; j<nstream; j++) flow_rate_sum += flow_rate[j];
    assert(flow_rate_sum!=0.);
    for (int j=0; j<nstream; j++) flow_rate[j] /= flow_rate_sum;
  }
  void time_step_pre(double time,int step) { }
  void time_step_post(double time,int step,
		      const vector<double> &gather_values) { 
    double sum_flow_rate_now=0., sum_flow_rate=0., p_avrg=0.;
    for (int j=0; j<nstream; j++) {
      flow_rate_now[j] = gather_values[gather_pos[j]];
      sum_flow_rate_now += flow_rate_now[j];
      p_avrg += pressure[j];
    }
    p_avrg /= double(nstream);
    PetscPrintf(PETSC_COMM_WORLD,
		"total flow rate: %g, avrg. press: %g\n",
		sum_flow_rate_now,p_avrg);
    for (int j=0; j<nstream; j++) {
      double new_pressure = pressure[j] + 
	flow_coef*(flow_rate_now[j]/sum_flow_rate_now-flow_rate[j]);
      PetscPrintf(PETSC_COMM_WORLD,
		  "stream %d, flow%d %g, desired f. %g, p%d %g, new_p%d %g\n",
		  j,flow_rate_now[j],j,flow_rate[j],j,pressure[j],j,new_pressure);
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
#endif
