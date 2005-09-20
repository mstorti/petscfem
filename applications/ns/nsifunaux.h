// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nsifunaux.h,v 1.1.2.1 2005/09/20 00:58:38 mstorti Exp $
#ifndef PETSCFEM_NSI_FUN_AUX_H
#define PETSCFEM_NSI_FUN_AUX_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** 
Some auxiliary functions commonly used by NSI applications
*/ 

// compute velocity dropping its projection in the axisymmetric direction  
void vel_axi(FastMat2 &u,FastMat2 &u_axi, int axi);

// compute h_supg for upwind 
double compute_h_supg(FastMat2 &u_axi,FastMat2 &dshapex, double velmod, double h_pspg);

// compute mixture density (rho_m)
double compute_rho_m(FastMat2 &rho_g_vp,FastMat2 &arho_g_vp,FastMat2 &alpha_g_vp, 
		     double &alpha_l, double rho_l, int nphases);

double compute_rho_m_2(FastMat2 &rho_g_vp,FastMat2 &arho_g_vp, FastMat2 &alpha_g_vp, 
		       double &alpha_l, double rho_l, int nphases, FastMat2 &grad_alpha_g_vp, 
		       FastMat2 &grad_rho_m,const int ndim);

double compute_rho_m_2_bcc(FastMat2 &rho_g_vp,FastMat2 &arho_g_vp, FastMat2 &alpha_g_vp, 
		       double &alpha_l, double rho_l, int nphases,const int ndim);

// compute disperse phase velocities
void compute_vel_g(const FastMat2 &v_m, FastMat2 &vslip_vp, const FastMat2 &vslip_user_vp, 
		   FastMat2 &vslip_m_vp, FastMat2 &v_g_vp, const FastMat2 &rho_g_vp, 
		   const double rho_m, const int g_dir, const FastMat2 &d_bubble_vp, 
		   const int nphases, const int use_modified_slag_vslip,
		   const FastMat2 &alpha_g_vp );

void compute_vel_g_2(const FastMat2 &V_v, FastMat2 &vslip_vp, const FastMat2 &vslip_user_vp, 
		     FastMat2 &vslip_m_vp, const FastMat2 &rho_g_vp,const double rho_m,
		     const double rho_l,const int g_dir,const FastMat2 &d_bubble_vp,
		     const int nphases,const int use_modified_slag_vslip,const double alpha_l,
		     const FastMat2 &alpha_g_vp,FastMat2 &v_g_vp,FastMat2 &grad_v_g_vp,
		     const FastMat2 &grad_rho_m,FastMat2 &grad_alpha_g_vp,const int ndim,
		     const FastMat2 &grad_V_v);

void compute_vel_g_2_bcc(const FastMat2 &V_v, FastMat2 &vslip_vp, const FastMat2 &vslip_user_vp, 
		     FastMat2 &vslip_m_vp, const FastMat2 &rho_g_vp,const double rho_m,
		     const double rho_l,const int g_dir,const FastMat2 &d_bubble_vp,
		     const int nphases,const int use_modified_slag_vslip,const double alpha_l,
			 const FastMat2 &alpha_g_vp,FastMat2 &v_g_vp,const int ndim);
#endif


