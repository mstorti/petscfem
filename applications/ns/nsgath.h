// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: nsgath.h,v 1.2.24.1 2004/08/22 23:01:11 mstorti Exp $
#ifndef NS_GATHERER_H
#define NS_GATHERER_H

extern int MY_RANK;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the force that the container is
    performing on the wall. 
    # force = - \int_\surf p * normal \dS #, where
    #force# is the computed force, #surf# is the
    surface of the elemset, #p# is pressure, #normal#
    is the unit vector \emph{exterior} to the fluid. The convention is
    that quad and triangular panels in 3D are numbered
    counterclock-wise when viewed from the exterior to the fluid. In
    2D linear elements are numbered following the countour curve in
    such a way as to leave the fluid at the left. 
    #\dS# is the diferential of surface. 
*/ 
class force_integrator : public gatherer {
private:
  /** Auxiliary vectors to store force, moment, the point about
      which to compute moments, the displacemente vector from the
      moment center to the actual integration point .
  */
  FastMat2 force,moment,x_center,dx;
  /// Flag to compute moments or not, number of dimensions
  int compute_moment, ndim_m;
public:
  /// perform several checks and initialization
  void init();
  /// set forces 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
  /// perform some cleaning
  void clean();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the flow rate at wall. 
    fixme:= agregar doc.
*/ 
class flow_rate_integrator : public gatherer {
private:
  FastMat2 Q;
  int ndim_m;
public:
  /// perform several checks and initialization
  void init();
  /// set forces 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the total height of the free surface (linear free surface
    b.c. version).
    fixme:= agregar doc.
*/ 
class free_surface_level_integrator : public gatherer {
private:
  int normal_dim;
public:
  /// perform several checks and initialization
  void init() { 
    assert(gather_length==1); 
    //o Position in gather vector
    int ierr;
    TGETOPTDEF_ND(thash,int,normal_dim,1);
    // if (!MY_RANK) printf("dim: %d\n",normal_dim);
  }
  /// set forces 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time) {
    pg_values[0] = wpgdet * xpg.get(normal_dim);
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the total volume occupied by the volume elemset. */ 
class volume_integrator : public gatherer {
private:
public:
  /// perform several checks and initialization
  void init() { assert(gather_length==1); }
  /// add volume
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time) {
    pg_values[0] = wpgdet;
  }
};

#endif
