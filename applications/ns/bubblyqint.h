// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: bubblyqint.h,v 1.2 2006/09/04 17:49:42 mstorti Exp $
#ifndef PETSCFEM_BUBBLYQINT_H
#define PETSCFEM_BUBBLYQINT_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class bubbly_flow_rate_integrator : public SurfGatherer {
private:
  FastMat2 vslip_user_vp, tmp, ;
  int nphases, ndim, g_dir;
public:
  void init();
  void set_ip_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &xpg,FastMat2 &n,double time);
  int vals_per_plane();
};

#endif
