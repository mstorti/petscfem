// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: advdfgth.h,v 1.1 2003/01/25 15:28:55 mstorti Exp $
#ifndef ADVDIF_GATHERER_H
#define ADVDIF_GATHERER_H

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

#endif
