// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: srfgath.h,v 1.4 2004/01/29 21:07:30 mstorti Exp $
#ifndef PETSCFEM_SRF_GATH_H
#define PETSCFEM_SRF_GATH_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Allows to compute integrals, or any kind of function
    that returns scalar values.
*/
class SurfGatherer : public Elemset { 
public:
  class SurfFunction {
  public:
    virtual void init(const TextHashTable *thash) { }
    virtual double f(const FastMat2 &x)=0;
    virtual void close() { }
    virtual ~SurfFunction()=0;
    static SurfFunction* factory(const TextHashTable *thash);
  };

private:
  SurfFunction *sf;

public: 
  SurfGatherer() : sf(NULL) { }
  ~SurfGatherer();

  int gather_length;

  /// The assemble function for the elemset. 
  ASSEMBLE_FUNCTION;
  /// The ask function for the elemset. 
  ASK_FUNCTION;

  /// Perform initialization
  void initialize();
  
  /** @name User defined callback functions
      @doc Deriving the class and defining these
      functions you can define several integrators.
  */
  //@{
  /// Called \textbf{before} the element loop. 
  virtual void init() { }

  /** Called for each element, before the Gauss point loop. 
      @param k (input) element number in the elemeset
  */ 
  virtual void element_hook(int k) {}

  /** Called at each Gauss point loop. Compute the contribution at
      each Gauss point. 
      @param pg_values (output) your contributions at the Gauss point
      @param u (input) the state of variables at time n+1
      @param uold (input) the state of variables at time n
      @param xpg (input) coordinates of the Gauss point
      @param n (input) normal to the surface (if #ndimel=ndim-1#)
      between real and master coordinates
      @param wpgdet (input) the weight of the Gauss point
      @param time (input) time value ($t^\alpha$)
      
  */ 
  virtual void set_ip_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &xpg,FastMat2 &n,double time)=0;

  /** Called \textbf{after} the element loop. May be used for
      clean-up operations. 
  */
  virtual void clean() {};
  //@}

  void handle_error(int error);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class field_surf_integrator : public SurfGatherer {
public:
  void set_ip_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &xpg,FastMat2 &n,double time);
};

#endif
