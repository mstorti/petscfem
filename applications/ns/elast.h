// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: elast.h,v 1.1.2.1 2001/10/29 14:34:41 mstorti Exp $

#ifndef ELASTICITY_H
#define ELASTICITY_H

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  elasticity : public adaptor { 
public: 
  double rho,E,nu;
  int ntens,nen;
  FastMat2 B,C,Jaco,iJaco,strain,stress,res_pg,mat_pg1,mat_pg2;
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/// 
class  elasticity_f : public adaptor { 
public: 
  double rho,E,nu,thickness,yield_strength;
  int ntens,nen;
  FastMat2 wpgdet_v,Jaco,iJaco,ehist_m,ehist_f,
    emass_m,gpcod_m,gpcod_f,epmtx_m,props_m,rmat1_m,
    shape_f,stran_m,stran_f,strsg_m,strsg_f,stra0_m,stra0_f,
    strs0_m,strs0_f,bmatx_m,bmatx_f,desig_m,dmatx_m,dmatx_f,
    dstra_m,presg_m,sgtot_m,sigma_m,tstra_m,tenod_m,dteno_m,
    pwoel_m,tgapl_m,preal_m,tenoi_m,xjacm_m,xjacm_f,
    bsbar_m,bsbar_f,xjacn_m,xjacn_f,fpchl_m,fpchl_f,
    hstif_m,hstif_f,
    dshapex_f; // To be improved and enhanced
  void init();
  void element_connector(const FastMat2 &xloc,
			 const FastMat2 &state_old,
			 const FastMat2 &state_new,
			 FastMat2 &res,FastMat2 &mat);
};

#endif
