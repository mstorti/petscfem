//__INSERT_LICENSE__
//$Id: tempfun2.cpp,v 1.4 2003/11/25 02:10:22 mstorti Exp $

#include <math.h>

#include "fem.h"
#include "getprop.h"
#include "dofmap.h"
#include "ampli.h"
#include "utils.h"
#include "util2.h"

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
    // ( #half_width = 0.83255 * sigma# )
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
      // Apply linear transformation to  time and amplitude
      double ttt = t_0 + t_scale * tt;
      double fff = A_0 + A_scale * ff;
      // Check that time vector is ordered. Compare with the last
      // entered value
      if (ntime && ttt <= t_v[ntime-1]) {
	PETSCFEM_ERROR("piecewise_linear: times must be a strictly increasing sequence\n"
		       "line: t=%f f=%f, previous time %f\n",tt,ff,(t_v[ntime-1]-t_0)/t_scale);
      }  
      // Load on vector
      t_v.push_back(ttt);
      f_v.push_back(fff);
      // increment counter
      ntime++;
    }
    fclose(fid);
    assert (ntime == t_v.size());
    assert (ntime>1);
    // Store pointers
    first = t_v.begin();
    end = t_v.end();
    t = &*t_v.begin();
    f = &*f_v.begin();
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
/** This time function represents
    #phi(t) = base + A * (4*(t-T0)*(T0+T-t)/T^2)^expo# if
    #T0<= t <= T0+T#, and #phi=0# otherwise.
 */
class smooth_impulse : public Amplitude {
private:
  double T0,T,A,expo,base;
  TextHashTable *thash;
public:
  void print() const {
    PetscPrintf(PETSC_COMM_WORLD,
		"\"smooth_impulse\" fixa amplitude, options table: \n");
    thash->print();
  }
  void init(TextHashTable *thash_) {
    thash = thash_;
    int ierr;
    //o The starting time
    TGETOPTDEF_ND(thash,double,T0,0.);
    //o The duration of the impulse
    TGETOPTDEF_ND(thash,double,T,1.);
    assert(T>0.);
    //o amplitud of impulse
    TGETOPTDEF_ND(thash,double,A,1.);
    //o exponent of function
    TGETOPTDEF_ND(thash,double,expo,2.);
    assert(expo>=0.);
    //o base value
    TGETOPTDEF_ND(thash,double,base,0.);
  }
  double eval(const TimeData *time_data) {
    double tt = double(* (const Time *) time_data);
    double xi;
    if (tt>=T0+T || tt<=T0) xi=0.;
    else xi = 4.*(tt-T0)*(T0+T-tt)/square(T);
    return base + A * pow(xi,expo);
  }
  ~smooth_impulse() { delete thash; }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class trigo_cos : public Amplitude {
protected:
  double mean_val, amplitude, omega, frequency, period, phase;
  TextHashTable *thash;
public:
  void print() const {
    PetscPrintf(PETSC_COMM_WORLD,
		"\"cos\" fixa amplitude, options table: \n");
    thash->print();
  }
  void init(TextHashTable *thash_) {
    thash = thash_;
    int ierr;
    SGETOPTDEF_ND(double,mean_val,0.);
    SGETOPTDEF_ND(double,amplitude,1.);
    SGETOPTDEF_ND(double,omega,0.);
    SGETOPTDEF_ND(double,frequency,0.);
    SGETOPTDEF_ND(double,period,0.);
    SGETOPTDEF_ND(double,phase,0.);
    if (omega==0. && frequency==0. && period==0) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "cos_function: not defined any of omega/frequency/period\n");
      PetscFinalize();
      exit(0);
    }
    if (frequency!=0.) omega = 2*M_PI*frequency;
    if (period!=0.) omega = 2*M_PI/period;
  }
  double eval(const TimeData *time_data) {
    double tt = double(* (const Time *) time_data);
    return mean_val + amplitude * cos(omega*tt+phase);
  }
  ~trigo_cos() { delete thash; }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class trigo_sin : public trigo_cos {
public:
  double eval(const TimeData *time_data) {
    double tt = double(* (const Time *) time_data);
    return mean_val + amplitude * sin(omega*tt+phase);
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Amplitude *Amplitude::factory(char *& label,
			      TextHashTable *t) {
  Amplitude *amp;
  if (!strcmp(label,"gaussian")) {
    amp = new gaussian;
  } else if (!strcmp(label,"piecewise_linear")) {
    amp = new piecewise_linear;
  } else if (!strcmp(label,"smooth_impulse")) {
    amp = new smooth_impulse;
  } else if (!strcmp(label,"cos")) {
    amp = new trigo_cos;
  } else if (!strcmp(label,"sin")) {
    amp = new trigo_sin;
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
