// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: ampli.h,v 1.5 2002/02/10 12:47:00 mstorti Exp $
#ifndef AMPLI_H
#define AMPLI_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This is the type of temporal functions. A TextHashTable stores
    various parameters, and the time (or time like data) is passed to
    the function. 
    @author M. Storti
*/ 
typedef double AmplitudeFunction(TextHashTable *thash,const TimeData *time_data,
				 void *& fun_data);

AmplitudeFunction smooth_ramp_function;
AmplitudeFunction ramp_function;
AmplitudeFunction sin_function;
AmplitudeFunction cos_function;
AmplitudeFunction piecewise_function;
AmplitudeFunction spline_function;
AmplitudeFunction spline_periodic_function;
typedef map<string,AmplitudeFunction *> FunctionTable;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** An amplitude is basically a pointer to a
    function that depends on some fixed parameters passed through
    a hash table and a variable parameter (typically time). 
    For instance function `sine'  with constant parameters  omega and
    phase, depending on time. 
    @author M. Storti
*/ 
class OldAmplitude : public Amplitude {
private:
  /// The key, i.e. the string identifying the function (e.g. `sine')
  char *amp_function_key;
  /** A table with parameters for the function (amplitude, starting
      time, for instance. )
  */
  TextHashTable *thash;
  /// A place where to store things
  void *fun_data;
public:
  /** A static table from the key to a pointer to the corresponding function. 
  */
  static FunctionTable *function_table;
  /// Constructor
  OldAmplitude(char *& s_,TextHashTable *tht_=NULL) :
    amp_function_key(s_), thash(tht_) {fun_data=NULL;};
  /// Eval the amplitude of the function at this time. 
  double eval(const TimeData *time_data);
  /// Adds an entry to the static table. 
  static void add_entry(const char * s,AmplitudeFunction *f);
  /// Initializes the function table. 
  static void initialize_function_table(void);
  /// Reads the table from a fstack
  void read_hash_table(FileStack *fstack);
  /// prints the amplitude entry. 
  void print(void) const; 
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class DLGeneric : public Amplitude {
private:
  TextHashTable *thash;
  typedef double EvalFun(double,void *);
  typedef void InitFun(TextHashTable *,void *&);
  typedef void ClearFun(void *);
  EvalFun *fun;
  InitFun *init_fun;
  ClearFun *clear_fun;
  void *handle;
  void *fun_data; // store data
public:
  DLGeneric() : fun_data(NULL) {}
  void print() const;
  void init(TextHashTable *thash_);
  virtual void clear() {};
  double eval(const TimeData *time_data);
  ~DLGeneric();
};

#endif
