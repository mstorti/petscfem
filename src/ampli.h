// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: ampli.h,v 1.11 2002/02/10 23:30:55 mstorti Exp $
#ifndef AMPLI_H
#define AMPLI_H

// Preprocessor macro `USE_DLEF' flags using dynamically loaded
// extended functions via `dlopen()' or not.
// If don't use this extensions then this file has little to give.

#ifdef USE_DLEF

#include <math.h>
#include <src/fem.h>
#include <src/getprop.h>

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
/// Generic amplitude function that dynamically loads functions
class DLGeneric : public Amplitude {
private:
  /// Options table
  TextHashTable *thash;
  /// Initialization function prototype
  typedef void InitFun(TextHashTable *,void *&);
  /// Evaluation function prototype
  typedef double EvalFun(double,void *);
  /// Cleanup function prototype
  typedef void ClearFun(void *);

  /// Initialization function
  InitFun *init_fun;
  /// Evaluation function
  EvalFun *eval_fun;
  /// Cleanup function
  ClearFun *clear_fun;

  /// Contains pointers for a given set of functions
  struct FunHandle {
    InitFun *init_fun;
    EvalFun *eval_fun;
    ClearFun *clear_fun;
  };

  /// This is the internal table that is kept for each file
  typedef map<string,FunHandle> FunTable;
  /** For each file we keep a handle (from `dlopen()') 
      table of pointers to functions
  */
  struct FileHandle {
    /// handle obtained from `dlopen()'
    void *handle;
    /// table of functions for each file
    FunTable *fun_table;
  };

  /** We keep a static table filename -> file_handle. 
      This is the generic type for this table.
  */
  typedef map<string,FileHandle> FileHandleTable;
  /// The actual `dlopen' handle for this instance
  void *handle;
  /// The actual handle for this instance
  FileHandle fh;
  /// The actual function handle for this instance
  FunHandle fuh;
  /** This generic pointer may be used to store internal values
      for the functions
  */
  void *fun_data;
  /// This is the actual table
  static FileHandleTable file_handle_table;

  string function_name,ext_filename;
public:
  /// Constructor (initializes `fun_data')
  DLGeneric() : fun_data(NULL) {}
  /// Prints information about the function
  void print() const;
  /// Initializes values 
  void init(TextHashTable *thash_);
  /// Clears the object 
  virtual void clear() {};
  /// Computes the value of the Dirichlet b.c. at the specified time
  double eval(const TimeData *time_data);
  /// Destructor
  ~DLGeneric();
};

// Usefule macros for defining extended functions

#define INIT_FUN extern "C" void init_fun(TextHashTable *thash,void *&fun_data)

#define INIT_FUN1(name) extern "C" \
        void name##_init_fun(TextHashTable *thash,void *&fun_data)

#define EVAL_FUN extern "C" double eval_fun(double t,void *fun_data)

#define EVAL_FUN1(name) extern "C" double name##_eval_fun(double t,void *fun_data)

#define CLEAR_FUN extern "C" void clear_fun(void *fun_data) 

#define CLEAR_FUN1(name) extern "C" void name##_clear_fun(void *fun_data) 

#define DEFINE_EXTENDED_AMPLITUDE_FUNCTION(fun_obj_class)	\
								\
INIT_FUN1(fun_obj_class) {					\
  fun_obj_class *fun_obj_ptr = new fun_obj_class;		\
  fun_data = fun_obj_ptr;					\
  fun_obj_ptr->init(thash);					\
}								\
								\
EVAL_FUN1(fun_obj_class) {					\
  fun_obj_class *fun_obj_ptr = (fun_obj_class *) fun_data;	\
  return fun_obj_ptr->eval_fun(t);				\
}								\
								\
CLEAR_FUN1(fun_obj_class) {					\
  fun_obj_class *fun_obj_ptr = (fun_obj_class *) fun_data;	\
  delete fun_obj_ptr;						\
}
#endif

#endif
