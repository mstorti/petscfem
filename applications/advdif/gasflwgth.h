// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: gasflwgth.h,v 1.1 2005/07/25 01:14:55 mstorti Exp $
#ifndef PETSCFEM_ADVDIF_GASFLWGTH_H
#define PETSCFEM_ADVDIF_GASFLWGTH_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the flow rate at wall. 
    fixme:= agregar doc.
*/ 
class gasflow_force_integrator : public gatherer {
private:
  FastMat2 F;
  int ndim_m;
public:
  /// perform several checks and initialization
  void init();
  /// set forces 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
};

#endif
