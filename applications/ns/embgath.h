// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: embgath.h,v 1.4 2002/08/07 11:47:26 mstorti Exp $
#ifndef EMBGATH_H
#define EMBGATH_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Surf2Vol : public GPdata {
private:
  int use_exterior_normal_m;
public:
  Surf2Vol(const char *geom,int ndim,int nel,
	   int npg,int mat_version=GP_NEWMAT,
	   int use_exterior_normal_a=0) 
    : GPdata(geom,ndim,nel,npg,mat_version),
  use_exterior_normal_m(use_exterior_normal_a) {}
  void map_mask(const int *map_fc,int *vicorow);
  int use_exterior_normal() { return use_exterior_normal_m; }
  virtual void face(int j,const int *&fc,const int *&vol)=0;
  virtual int nfaces(int &nel_surf,int &nel_vol)=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Quad2Hexa : public Surf2Vol {
private:
  static const int faces[][8];
  int vol[8], vol_r[8];
  int this_face[4];
public:
  Quad2Hexa(const char *geom,int ndim,int nel,
	   int npg,int mat_version=GP_NEWMAT,
	    int use_exterior_normal_m=0) 
    : Surf2Vol(geom,ndim,nel,npg,mat_version,
	       use_exterior_normal_m) { }
  void face(int j,const int *&fc,const int *&vol);
  int nfaces(int &nel_surf,int &nel_vol) { nel_surf=4; nel_vol=8; return 24; }
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Allows to compute integrals, or any kind of function
    that returns scalar values.
*/
class embedded_gatherer : public Elemset { 
public: 
  void initialize();
  int gather_length;
  /// This should not be defined by the user
  ASSEMBLE_FUNCTION;
  /// Return which "jobinfos" will process
  ASK_FUNCTION;
  
  /** @name User defined callback functions
      @doc Deriving the class and defining these
      functions you can define several integrators.
  */
  //@{
  /** Number of surface nodes, desired number of
      nodes in the volume element. */
  virtual void surface_nodes(int &nel_surf,int &nel_vol)=0;
  /// Called \textbf{before} the element loop. 
  virtual void init() { }

  /** Called for each element, \emph{before} the Gauss point loop. 
      @param k (input) element number in the elemeset
  */ 
  virtual void element_hook(int k) {}

  /** Called at each Gauss point loop. Compute the contribution at
      each Gauss point. 
      @param pg_values (output) your contributions at the Gauss point
      @param u (input) the field state at time n+1
      @param uold (input) the field state at time n
      @param grad_u (input) the gradient of field state of variables at time n+1
      @param grad_uold (input) the gradient of field state of variables at time n
      @param xpg (input) coordinates of the Gauss point
      @param n (input) normal to the surface (if #ndimel=ndim-1#)
      between real and master coordinates
      @param wpgdet (input) the weight of the Gauss point
      @param time (input) time value ($t^\alpha$)
  */ 
  virtual void set_pg_values(vector<double> &pg_values,FastMat2 &u,
			     FastMat2 &uold,FastMat2 &grad_u, FastMat2 &grad_uold, 
			     FastMat2 &xpg,FastMat2 &n,
			     double wpgdet,double time)=0;

  /** Called \textbf{after} the element loop. May be used for
      clean-up operations. 
  */
  virtual void clean() {};
  //@}
};

class visc_force_integrator : public embedded_gatherer { 
public:
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &grad_u, FastMat2 &grad_uold, 
		     FastMat2 &xpg,FastMat2 &n,
		     double wpgdet,double time);
  void surface_nodes(int &nel_surf,int &nel_vol);
};

#endif
