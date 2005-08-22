// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: advdfgth.h,v 1.3 2005/08/22 19:42:41 mstorti Exp $
#ifndef PETSCFEM_ADVDIF_GATHERER_H
#define PETSCFEM_ADVDIF_GATHERER_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes the flow rate at wall. 
    fixme:= agregar doc.
*/ 
class flow_rate_integrator : public gatherer {
private:
  FastMat2 Q;
  int ndim_m, use_mass_rate;
public:
  /// perform several checks and initialization
  void init();
  /// set forces 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		     double wpgdet,double time);
};

#endif
