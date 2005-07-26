// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: gasflwgth.h,v 1.2 2005/07/26 15:36:04 mstorti Exp $
#ifndef PETSCFEM_ADVDIF_GASFLWGTH_H
#define PETSCFEM_ADVDIF_GASFLWGTH_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the flow rate at wall. 
    fixme:= agregar doc.
*/ 
class gasflow_force_integrator : public gatherer {
private:
  FastMat2 F,x0,M,dx;
  int ndim_m, comp_moments, dim_mom;
public:
  /// perform several checks and initialization
  void init();
  /// set forces 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
};

#endif
