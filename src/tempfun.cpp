//__INSERT_LICENSE__
//$Id: tempfun.cpp,v 1.11 2002/04/10 19:25:26 mstorti Exp $

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
  add_entry("sin",&sin_function);
  add_entry("cos",&cos_function);
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
  } else if (end_time>start_value && time>end_time) {
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class gaussian : public Amplitude {
private:
  double base, A, sigma, t0;
  TextHashTable *thash;
public:
  void print() const {
    PetscPrintf(PETSC_COMM_WORLD,
		"\"gaussian\" fixa amplitude, options table: \n");
    thash->print();
  }
  void init(TextHashTable *thash_) {
    thash = thash_;
    int ierr;
    double peak, half_width;

    //o The base value for the Gaussian
    TGETOPTDEF_ND(thash,double,base,0.);
    assert(ierr==0);
    
    //o The absolute peak value of the Gaussian
    TGETOPTDEF_ND(thash,double,peak,base+1.);
    assert(ierr==0);

    //o The height (with reference to the base) of the Gaussian
    TGETOPTDEF_ND(thash,double,A,peak-base);
    assert(ierr==0);

    //o The width of the Gaussian, i.e. the time value such that the
    // relative amplitude falls from the peak value to one half
    // (#half_width = 0.83255 * sigma#)
    TGETOPTDEF_ND(thash,double,half_width,0.5);
    assert(ierr==0);

    //o The standard deviation of the Gaussian
    TGETOPTDEF_ND(thash,double,sigma,half_width/0.83255);
    assert(ierr==0);

    //o The time instant at which the peak of the Gaussian occurs 
    TGETOPTDEF_ND(thash,double,t0,0.);
    assert(ierr==0);
  }
  double eval(const TimeData *time_data) {
    double t = double(* (const Time *) time_data);
    double f = base + A * exp(-square((t-t0)/sigma));
    return f;
  }
  ~gaussian() { delete thash; }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class piecewise_linear : public Amplitude {
private:
  int ntime, periodic;
  vector<double> t_v,f_v;
  vector<double>::iterator first,end;
  double *t,*f, tini,tend;
  TextHashTable *thash;
public:
  void print() const {
    PetscPrintf(PETSC_COMM_WORLD,
		"\"piecewise_linear\" fixa amplitude, options table: \n");
    thash->print();
    PetscPrintf(PETSC_COMM_WORLD,"number of intervals: %f\n",ntime);
  }
  void init(TextHashTable *thash_) {
    double tt,ff;
    thash = thash_;
    int ierr;

    //o The name of the file
    TGETOPTDEF_S(thash,string,filename,piecewise_linear_data.dat);
    assert(ierr==0);

    //o Shift in the time linear transformation.
    // $t = t_0 + t_scale * tau$. $t$ is time and $\tau$
    // is the pseudo-time entered in the second column of
    // the file. 
    TGETOPTDEF(thash,double,t_0,0.);

    //o Temporal scale in the time linear transformation.
    TGETOPTDEF(thash,double,t_scale,1.);

    //o Shift in the amplitud value linear transformation.
    // $A = A_0 + A_scale * \alpha$, where $A$ is the amplitude to
    // apply to the diplacements and $\alpha$ are the values entered
    // in the file. 
    TGETOPTDEF(thash,double,A_0,0.);

    //o Temporal scale in the time linear transformation.
    TGETOPTDEF(thash,double,A_scale,1.);

    // Read a table of t/f from file
    // fixme:= read filename from thash
    // string filename = "pcwzlin.dat.tmp";
    FILE *fid = fopen(filename.c_str(),"r");
    assert(fid);
    ntime=0;
    while (1) {
      int nread = fscanf(fid,"%lf %lf",&tt,&ff);
      if (nread==EOF) break;
      assert(nread==2);
      // increment counter
      ntime++;
      // Apply linear transformation to  time and amplitude
      tt = t_0 + t_scale*(tt-t_0);
      ff = A_0 + A_scale*(ff-A_0);
      // Check that time vector is ordered. Compare with the last
      // entered value
      if (ntime>=2) assert(tt > t_v[ntime-2]);
      // Load on vector
      t_v.push_back(tt);
      f_v.push_back(ff);
    }
    fclose(fid);
    assert (ntime == t_v.size());
    assert (ntime>1);
    // Store pointers
    first = t_v.begin();
    end = t_v.end();
    t = t_v.begin();
    f = f_v.begin();
    tini = t_v[0];
    tend = t_v[ntime-1];
  }
  double eval(const TimeData *time_data) {
    vector<double>::iterator k;
    double tt = double(* (const Time *) time_data);
    assert(tt>=tini);
    assert(tt<=tend);
    int j1=0, j2=ntime-2, j;
    while (1) {
      if (j1==j2) {
	break;
      } else if (j2==j1+1) {
	if (tt>t[j2]) j1=j2;
	break;
      }
      // printf("j1 %d, t: %f,  j2: %d, t=%f\n",
      // j1,t[j1],j2,t[j2]);
      j=(j1+j2+1)/2;
      if (tt>t[j]) {
	j1=j;
      } else {
	j2=j;
      }
      // printf("j, %d, t: %f\n",j,t[j]);
    }
    double slope = (f[j1+1]-f[j1])/(t[j1+1]-t[j1]);
    return f[j1]+slope*(tt-t[j1]);
      
  }
  ~piecewise_linear() { 
    delete thash; 
    // THis is not needed, but anyway...
    t_v.clear();
    f_v.clear();
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Amplitude *Amplitude::factory(char *& label,
			      TextHashTable *t=NULL) {
  Amplitude *amp;
  if (!strcmp(label,"gaussian")) {
    amp = new gaussian;
  } else if (!strcmp(label,"piecewise_linear")) {
    amp = new piecewise_linear;
  } else if (!strcmp(label,"dl_generic")) {
#ifdef USE_DLEF
    amp = new DLGeneric;
#else
    PetscPrintf(PETSC_COMM_WORLD,
		"This version is not compiled with dynamically"
		" loaded extended functions!!\n"
		" Enable the 'USE_DLEF' flag and recompile.\n");
    assert(0);
#endif
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
		"Not known fixa_amplitude \"%s\"\n",label);
    assert(0);
  }
  amp->init(t);
  delete t; // If it wasn't deleted by `init'
  return amp;
}
  
