//__INSERT_LICENSE__
//$Id: tempfun.cpp,v 1.15 2003/07/13 17:19:53 mstorti Exp $

#include <math.h>

#include "fem.h"
#include "getprop.h"
#include "dofmap.h"
#include "ampli.h"
#include "utils.h"
#include "util2.h"

using namespace std;

FunctionTable *OldAmplitude::function_table=NULL;

// Useful for definitions inside functions defining a `fun-data' class
// object 
#define FGETOPTDEF(type,name,default) \
        sd->name = default; \
        ierr = get_##type(thash,#name,&sd->name,1); \
        PFEMERRCA(ierr,"Error getting option \"" #name "\"\n") 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void OldAmplitude::initialize_function_table(void) {
  if (function_table) return; // already loaded may be
  function_table = new FunctionTable;
  add_entry("smooth_ramp",&smooth_ramp_function);
  add_entry("ramp",&ramp_function);
  // add_entry("sin",&sin_function);
  // add_entry("cos",&cos_function);
  add_entry("piecewise",&piecewise_function);
  add_entry("spline",&spline_function);
  add_entry("spline_periodic",&spline_periodic_function);
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double cos_function(TextHashTable *thash,const TimeData *t,
		    void *& fun_data) {
  int ierr;

  SGETOPTDEF(double,mean_val,0.);
  SGETOPTDEF(double,amplitude,1.);
  SGETOPTDEF(double,omega,-1.);
  SGETOPTDEF(double,frequency,-1.);
  SGETOPTDEF(double,period,-1.);
  SGETOPTDEF(double,phase,0.);

  if (omega<0. && frequency<0. && period<0) {
    PetscPrintf(PETSC_COMM_WORLD,
		"cos_function: not defined any of omega/frequency/period\n");
    PetscFinalize();
    exit(0);
  }
  if (frequency>0.) omega = 2*M_PI*frequency;
  if (period>0.) omega = 2*M_PI/period;
//   double time = time->time();
  double time = *(Time *)t;
  return mean_val+amplitude*cos(omega*time+phase);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double sin_function(TextHashTable *thash,const TimeData *t,
		    void *& fun_data) {
  int ierr;

  SGETOPTDEF(double,mean_val,0.);
  SGETOPTDEF(double,amplitude,1.);
  SGETOPTDEF(double,omega,-1.);
  SGETOPTDEF(double,frequency,-1.);
  SGETOPTDEF(double,period,-1.);
  SGETOPTDEF(double,phase,0.);

  if (omega<0. && frequency<0. && period<0) {
    PetscPrintf(PETSC_COMM_WORLD,
		"sin_function: not defined any of omega/frequency/period\n");
    PetscFinalize();
    exit(0);
  }
  if (frequency>0.) omega = 2*M_PI*frequency;
  if (period>0.) omega = 2*M_PI/period;
//   double time = time_data->time();
  double time = *(Time *)t;
  //  double time=*t;
  double val=mean_val+amplitude*sin(omega*time+phase);
  // printf("time %f, amp: %f\n",time,val);
  return val;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double smooth_ramp_function(TextHashTable *thash,const TimeData *t,
		    void *& fun_data) {
  int ierr;
  SGETOPTDEF(double,switch_time,0.);
  SGETOPTDEF(double,time_scale,-1.);
  SGETOPTDEF(double,start_value,0.);
  SGETOPTDEF(double,end_value,1.);
  if (time_scale<0.) {
    PetscPrintf(PETSC_COMM_WORLD,
		"smooth_ramp_function: You mut specify '"
		"a positive `time_scale' parameter");
    PetscFinalize();
    exit(0);
  }
//   double time = time_data->time();
  double mean_val = (start_value+end_value)/2.;
  double amplitude = (end_value-start_value)/2.;
//   double time=*t;
  double time = *(Time *)t;
  double val = mean_val + amplitude *
    tanh((time-switch_time)/time_scale);
  return val;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double ramp_function(TextHashTable *thash,const TimeData *t,
		    void *& fun_data) {
  int ierr;
  SGETOPTDEF(double,start_time,0.);
  SGETOPTDEF(double,start_value,0.);
  SGETOPTDEF(double,slope,0.);
  SGETOPTDEF(double,end_time,start_time);
  SGETOPTDEF(double,end_value,start_value);
  if (slope==0.)
    slope = (end_value-start_value)/(end_time-start_time);

//   double time = *(double *)time_data;
  double time = *(Time *)t;
//   double time=*t;
  if (time<start_time) {
    return start_value;
  } else if (end_time>start_time && time>end_time) {
    return end_value;
  } else {
    return start_value+slope*(time-start_time);
  }
}

// In this structure we store 
struct piecewise_data {
  double *time_vals,*ampl_vals;
  int npoints;
  double final_time;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double piecewise_function(TextHashTable *thash,const TimeData *t,
			  void *& fun_data) {
  int ierr,npoints;
  double *time_vals,*ampl_vals,final_time;
  piecewise_data *pd;

  if (fun_data==NULL) {
    pd = new piecewise_data;
    fun_data = (void *)pd;

    SGETOPTDEF_ND(int,npoints,2);
    assert(npoints>=2);
    pd->npoints = npoints;

    pd->time_vals = new double[npoints];
    pd->ampl_vals = new double[npoints];

    ierr = get_double(thash,"time_vals",pd->time_vals,0,npoints); CHKERRA(ierr);
    ierr = get_double(thash,"ampl_vals",pd->ampl_vals,0,npoints); CHKERRA(ierr);

    SGETOPTDEF_ND(double,final_time,pd->time_vals[0]);
    pd->final_time = final_time;
  }
  pd = (piecewise_data *)fun_data;
  npoints = pd->npoints;
  time_vals = pd->time_vals;
  ampl_vals = pd->ampl_vals;
  final_time = pd->final_time;

  double time = *(Time *)t;
  double null_value = 0.;
  double start_time = time_vals[0];
  double end_time = time_vals[npoints-1];

  if (time<start_time) {
    return null_value;
  } else if (final_time>start_time && time>final_time) {
    return null_value;
  } else {
    // double tstar = start_time + (time-start_time) % (end_time-start_time);
    double tstar = start_time + fmod(time-start_time,end_time-start_time);
    double t1,t2,f1,f2;
    for (int k=0; k<npoints-1; k++) {
      t1 = time_vals[k];
      t2 = time_vals[k+1];
      f1 = ampl_vals[k];
      f2 = ampl_vals[k+1];
      if (tstar<=t2 && tstar>=t1) break;
    }
    if (t1==t2) {
      return 0.5*(f1+f2);
    } else {
      return f1 + (tstar-t1)/(t2-t1)*(f2-f1);
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// In this structure we store 
struct spline_data {
  double *time_vals,*ampl_vals;
  double final_time;
  int npoints,periodic;
  double *b,*c,*d; // spline coefficients
  double *bi,*ci,*di; // spline coefficients for the odd part
  
};

// Routines downloaded from Netlib
extern "C" {
  int spline_(int *n, double *x, double *y, double *b, double *c, double *d);
  double seval_(int *n, double *u, double *x,
	       double *y, double *b, double *c, double *d);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "AmplitudeFunction spline_function"
double spline_function(TextHashTable *thash,const TimeData *t,
		       void *& fun_data) {

  int ierr,npoints,periodic;
  double *time_vals,*ampl_vals,final_time;
  spline_data *sd;

  if (fun_data==NULL) {
    sd = new spline_data;
    fun_data = (void *)sd;

    SGETOPTDEF_ND(int,npoints,2);
    assert(npoints>=2);
    sd->npoints = npoints;

    sd->time_vals = new double[npoints];
    sd->ampl_vals = new double[npoints];

    sd->b = new double[npoints];
    sd->c = new double[npoints];
    sd->d = new double[npoints];

    ierr = get_double(thash,"time_vals",sd->time_vals,0,npoints); CHKERRA(ierr);
    ierr = get_double(thash,"ampl_vals",sd->ampl_vals,0,npoints); CHKERRA(ierr);

    SGETOPTDEF_ND(double,final_time,sd->time_vals[0]);
    SGETOPTDEF_ND(int,periodic,0);
    sd->final_time = final_time;

    spline_(&npoints,sd->time_vals,sd->ampl_vals,sd->b,sd->c,sd->d);
  }

  sd = (spline_data *)fun_data;
  npoints = sd->npoints;
  time_vals = sd->time_vals;
  ampl_vals = sd->ampl_vals;
  final_time = sd->final_time;

  double time = *(Time *)t;
  double null_value = 0.;
  double start_time = time_vals[0];
  double end_time = time_vals[npoints-1];

  if (!periodic) 
    final_time = time_vals[npoints-1];

  if (time<start_time) {
    return null_value;
  } else if (final_time>start_time && time>final_time) {
    return null_value;
  } else {
    // double tstar = start_time + (time-start_time) % (end_time-start_time);
    double tstar = start_time + fmod(time-start_time,end_time-start_time);
    double fspl = seval_(&npoints, &tstar, time_vals, ampl_vals, sd->b, sd->c, sd->d);
    return fspl;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// In this structure we store 
struct spline_periodic_data {
  double start_time,period;
  int npoints;
  double *feven,*fodd,*fun;    // function values and odd/even decomp.
  double *beven,*ceven,*deven; // spline coefficients
  double *bodd,*codd,*dodd;    // spline coefficients for the odd part
  double *x;                   // this is redundant but needed by the
			       // fortran spline functions
  double eval(double time);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double spline_periodic_function(TextHashTable *thash,const TimeData *t,
		       void *& fun_data) {

  int ierr;
  spline_periodic_data *sd;

  if (fun_data==NULL) {
    sd = new spline_periodic_data;
    fun_data = (void *)sd;

    FGETOPTDEF(double,period,0);
    FGETOPTDEF(int,npoints,2);
    FGETOPTDEF(double,start_time,0.);
    int npoints=sd->npoints;
    assert(npoints>=2);
    assert(npoints % 2==1);
    int nphalf = npoints/2+1;

    sd->beven = new double[nphalf];
    sd->ceven = new double[nphalf];
    sd->deven = new double[nphalf];

    sd->bodd = new double[nphalf];
    sd->codd = new double[nphalf];
    sd->dodd = new double[nphalf];

    sd->fun = new double[npoints];
    sd->feven = new double[nphalf];
    sd->fodd = new double[nphalf];
    sd->x = new double[nphalf];

    ierr = get_double(thash,"ampl_vals",sd->fun,0,npoints);
    CHKERRA(ierr);

    // For periodic problems convert time to phase and extend
    // symmetric 
    double delta_phi = M_PI/double(nphalf-1);
    sd->fodd[0] = (sd->fun[1] - sd->fun[npoints-2])/(2.0*sin(delta_phi));
    sd->fodd[nphalf-1] = (sd->fun[nphalf-2] - sd->fun[nphalf])/(2.0*sin(delta_phi));
    for (int j=0; j<nphalf; j++) {
      double phase = double(j)*M_PI/double(nphalf-1);
      sd->x[j] = (1.0-cos(phase))/2.0;
      sd->feven[j] = 0.5*(sd->fun[j]+sd->fun[npoints-j-1]);
      if (j!=0 && j!=nphalf-1) 	
	sd->fodd[j] =
	  0.5*(sd->fun[j]-sd->fun[npoints-j-1])/sin(phase);
    }

    spline_(&nphalf,sd->x,sd->feven,sd->beven,sd->ceven,sd->deven);
    spline_(&nphalf,sd->x,sd->fodd,sd->bodd,sd->codd,sd->dodd);

  }

  sd = (spline_periodic_data *)fun_data;
  double time = *(Time *)t;
  double ampl=sd->eval(time);
  return ampl;

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double spline_periodic_data::eval(double time) {

  double phase = fmod(time-start_time,period)/period;
  if (phase<0) phase += 1;
  phase *= 2.0*M_PI;
  double xx = (1.0-cos(phase))/2.0;
  int nphalf = npoints/2+1;
  double fe = seval_(&nphalf, &xx, x, feven, beven, ceven, deven);
  double fo = seval_(&nphalf, &xx, x, fodd, bodd, codd, dodd);
  double fspl = fe + fo*sin(phase);

#if 0 // debug
  static double old_time=-1;
  if (time>old_time) {
    printf("time,xx,fe,fo,phase,fspl: %f %f %f %f %f %f\n",
	   time,xx,fe,fo,phase,fspl);
    old_time = time;
  }
#endif 

  return fspl;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void OldAmplitude::read_hash_table(FileStack *fstack) {
  ::read_hash_table(fstack,thash);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void OldAmplitude::add_entry(const char * s,
			  AmplitudeFunction *f) {
  pair<string,AmplitudeFunction *> pp(s,f);
  function_table->insert(pp);
  // function_table.insert(pair(string(s),f));
  // function_table[string(s)] = f;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void OldAmplitude::print(void) const {
  thash->print();
  PetscPrintf(PETSC_COMM_WORLD,"amp_function_key: %s\n",amp_function_key);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double OldAmplitude::eval(const TimeData *time_data)  {
  assert(function_table!=NULL);
  AmplitudeFunction *func;
  FunctionTable::iterator it;
  it = function_table->find(amp_function_key);
  if (it==function_table->end()) {
    PetscPrintf(PETSC_COMM_WORLD,
		"Not found key \"%s\" in function table\n",
		amp_function_key);
    PetscFinalize();
    exit(0);
  }
  func = it->second;
  double amplitude;
  amplitude = (*func)(thash,time_data,fun_data);
  return amplitude;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Amplitude * Amplitude::old_factory(char *& label,FileStack &fstack) {
  FunctionTable::iterator it;
  if (!OldAmplitude::function_table) OldAmplitude::initialize_function_table();
  it = OldAmplitude::function_table->find(label);
  if (it == OldAmplitude::function_table->end()) {
    return NULL;
  } else {
    OldAmplitude *amp = new OldAmplitude(label);
    amp->read_hash_table(&fstack);
    return amp;
  }
} 

