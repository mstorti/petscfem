//__INSERT_LICENSE__
//$Id: gasflow.cpp,v 1.2 2003/01/26 15:06:36 mstorti Exp $
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
  double flow_resistance,char_time;
  string outlet_name;
public:
  double pressure;
  void init(Mesh &mesh,Dofmap &dofmap,
	    TextHashTableFilter *o,const char *name) { 
    TextHashTable *options = (TextHashTable *)TextHashTable::find(name);
    int ierr;
    TGETOPTDEF_ND(options,double,flow_resistance,0.);
    TGETOPTDEF_ND(options,double,char_time,0.);
    assert(flow_resistance>0.);
    assert(char_time>0.);
    TGETOPTDEF(options,double,initial_pressure,0.);
    assert(initial_pressure>0.);
    pressure = initial_pressure;
    outlet_name = string(name);
    printf("registering gasflow_hook: %s\n",outlet_name.c_str());
    assert(gasflow_hook_table.find(outlet_name)
	   ==gasflow_hook_table.end());
    gasflow_hook_table[outlet_name] = this;
  }
  void time_step_pre(double time,int step) { }
  void time_step_post(double time,int step,
		      const vector<double> &gather_values) { 
  }
  void close() { }
};

DL_GENERIC_HOOK(gasflow_hook);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class flow_controller : public DLGenericTmpl {
private:
  gasflow_hook *my_gasflow_hook;
public:
  flow_controller() : my_gasflow_hook(NULL) {}
  void init(TextHashTable *thash) { 
    int ierr;
    TGETOPTDEF_S(thash,string,outlet_name,<none>);
    printf("Looking for hook %s\n",outlet_name.c_str());
    gasflow_hook_table_t::iterator q = 
      gasflow_hook_table.find(outlet_name);
    assert(q!=gasflow_hook_table.end());
    my_gasflow_hook = q->second;
  }
  double eval(double) { return my_gasflow_hook->pressure;  }
};

DEFINE_EXTENDED_AMPLITUDE_FUNCTION2(flow_controller);
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
