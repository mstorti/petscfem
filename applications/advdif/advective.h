// -*- mode: c++ -*-
 
/*
  This file belongs to he PETSc - FEM package a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/
 
//#define CHECK_JAC // Computes also the FD Jacobian for debugging
 
#ifndef ADVECTIVE_H
#define ADVECTIVE_H

/** The jacobians of the flux functions. This is an array
    of ndim matrices of ndof x ndof entries each. 
*/
#define AJAC(jd) (*A_jac[(jd)-1])

/// Standard arguments for the typical flux function. Newmat version. 
#define FLUX_FUN_ARGS const RowVector &U,int ndim,const Matrix &iJaco,		\
	      Matrix &H, Matrix &grad_H,					\
	      Matrix &flux,vector<Matrix *> A_jac,				\
	      Matrix &A_grad_U, Matrix &grad_U,					\
	      Matrix &G_source,Matrix &tau_supg, double &delta_sc,		\
	      double &lam_max,							\
	      TextHashTable *thash,						\
              Matrix &nor,Matrix &lambda,Matrix &Vr,Matrix &Vr_inv,		\
              double *propel,							\
	      void *user_data,int options, int &start_chunk, int & ret_options

#define FLUX_FUN_ARGS_GENER(ML_) const ML_ &U,			\
              int ndim,const ML_ &iJaco,			\
	      ML_ &H, ML_ &grad_H,				\
	      ML_ &flux,vector<ML_ *> A_jac,			\
	      ML_ &A_grad_U, ML_ &grad_U,			\
	      ML_ &G_source,ML_ &tau_supg, double &delta_sc,	\
	      double &lam_max,					\
	      TextHashTable *thash,				\
              ML_ &nor,ML_ &lambda,ML_ &Vr,ML_ &Vr_inv,		\
              double *propel,					\
	      void *user_data,int options,			\
              int &start_chunk, int & ret_options

#define AD_FLUX_FUN_ARGS const  FastMat2  &U,					\
	      int ndim, const  FastMat2  &iJaco,  FastMat2  &H,			\
	      FastMat2  &grad_H,  FastMat2  &flux, FastMat2 &fluxd,		\
	      FastMat2  &A_jac,  FastMat2  &A_grad_U,				\
	      FastMat2  &grad_U,  FastMat2  &G_source,				\
              FastMat2  &D_jac, FastMat2 &C_jac,				\
	      FastMat2  &tau_supg, double &delta_sc, double &lam_max,		\
	      TextHashTable *thash,  FastMat2  &nor, FastMat2  &lambda,		\
	      FastMat2  &Vr, FastMat2  &Vr_inv, double *propel,			\
	      void *user_data,int options,int &start_chunk, int & ret_options

#define FLUX_FUN_ARGS_FM FLUX_FUN_ARGS_GENER(FastMat)

/// Standard arguments for the typical flux function (in the function call).
#define FLUX_FUN_CALL_ARGS_GENER(ML_) U##ML_,ndim,				\
              iJaco##ML_,H##ML_,grad_H##ML_,					\
              flux##ML_,A_jac##ML_,A_grad_U##ML_,				\
	      grad_U##ML_,G_source##ML_,					\
              tau_supg##ML_,delta_sc,lam_max,thash,				\
              nor##ML_,lambda##ML_,Vr##ML_,Vr_inv##ML_,propel,user_data,	\
	      options,start_chunk,ret_options

#define FLUX_FUN_CALL_ARGS_FM FLUX_FUN_CALL_ARGS_GENER(_FM)

#define FLUX_FUN_CALL_ARGS_FM2 FLUX_FUN_CALL_ARGS_GENER(_FM2)

#define FLUX_FUN_CALL_ARGS U,ndim,iJaco,H,grad_H,flux,A_jac,A_grad_U,	\
	      grad_U,G_source,tau_supg,delta_sc,lam_max,thash,		\
              nor,lambda,Vr,Vr_inv,propel,user_data,			\
	      options,start_chunk,ret_options
	      
/// Type of flux function. 
typedef int AdvDifFluxFunction (AD_FLUX_FUN_ARGS);

#define ADVDIFFF_ARGS const NewElemset *elemset,const  FastMat2  &U,		\
	   int ndim, const  FastMat2  &iJaco,  FastMat2  &H,			\
	   FastMat2  &grad_H,  FastMat2  &flux, FastMat2 &fluxd,		\
	   FastMat2  &A_jac,  FastMat2  &A_grad_U,				\
	   FastMat2  &grad_U,  FastMat2  &G_source,				\
	   FastMat2  &D_jac,	FastMat2  &C_jac,				\
	   FastMat2  &tau_supg, double &delta_sc, double &lam_max,		\
	   FastMat2  &nor, FastMat2  &lambda,					\
	   FastMat2  &Vr, FastMat2  &Vr_inv, double *propel,			\
	   void *user_data,int options,int &start_chunk, int & ret_options

/// Generic FastMat2 matrix
#if 0
class FastMat2Shell {
public:
  /// Performs A = (*this) * B;
  virtual apply(FastMat2 & A,FastMat2 & B) const;
}
#endif
typedef void FastMat2Shell(FastMat2 & A,FastMat2 & B);

#define COMPUTE_FLUX_ARGS							\
		const FastMat2 &U,						\
		const FastMat2 &iJaco, FastMat2 &H,				\
		FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,		\
		FastMat2 &A_grad_U,						\
		FastMat2 &grad_U, FastMat2 &G_source,				\
		FastMat2 &D_jac, FastMat2 &C_jac,				\
		FastMat2 &tau_supg, double &delta_sc,				\
		double &lam_max,						\
		FastMat2 &nor, FastMat2 &lambda,				\
		FastMat2 &Vr, FastMat2 &Vr_inv,					\
		int options

// This is the flux function for a given physical problem. 
class AdvDifFF {
private:
  // The list of variables to be logarithmically transformed
  vector<int> log_vars_v;
public:
  virtual int operator() (ADVDIFFF_ARGS) {
    PetscPrintf(PETSC_COMM_WORLD,"Undefined flux function\n");
    return 0;
  }
  virtual void get_log_vars(const NewElemset *elemset,int &nlog_vars, 
			    const int *& log_vars);
};

class AJac {
public:
  virtual FastMat2Shell comp_A_grad_U=0,comp_A_grad_N=0;
  virtual void comp_flux(FastMat2 & A,FastMat2 & B) =0 ;
};

class NewAdvDif;

// This is the flux function for a given physical problem. 
class NewAdvDifFF {
private:
  // The list of variables to be logarithmically transformed
  vector<int> log_vars_v;
public:
  const NewAdvDif *elemset;
  AJac *a_jac;
  NewAdvDifFF(NewAdvDif *elemset_) : elemset(elemset_) {};
  virtual void start_chunk(int ret_options) =0;
  virtual void element_hook(ElementIterator &element) =0;
  virtual void compute_flux(COMPUTE_FLUX_ARGS) =0;
  virtual void get_log_vars(int &nlog_vars,const int *& log_vars);
};

/** The class AdvDif is a NewElemset class plus a
    advdif flux function object.
*/
class AdvDif : public NewElemset { 
public:
  NewAssembleFunction new_assemble;
  ASK_FUNCTION;
  AdvDifFF *adv_diff_ff;
};

class BcconvAdv : public NewElemset { 
public:
  NewAssembleFunction new_assemble;
  ASK_FUNCTION;
  AdvDifFF *adv_diff_ff;
};

enum flux_fun_opt {
  DEFAULT     = 0x0000,
  COMP_SOURCE = 0x0001,
  COMP_UPWIND = 0x0002,
  COMP_EIGENV = 0x0004,
  SCALAR_TAU  = 0x0008,
};

/** The class AdvDif is a NewElemset class plus a
    advdif flux function object.
*/
class NewAdvDif : public NewElemset { 
  NewAdvDifFF *adv_diff_ff;
public:
  int ndim,ndof;
  NewAdvDif(NewAdvDifFF *adv_diff_ff_) : adv_diff_ff(adv_diff_ff_) {};
  NewAssembleFunction new_assemble;
  ASK_FUNCTION;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define ADVECTIVE_ELEMSET(name)				\
FluxFunction flux_fun_##name;				\
							\
class volume_##name : public Advective {		\
public:							\
  volume_##name() {flux_fun=&flux_fun_##name;};		\
};							\
							\
class absorb_##name : public Absorb {			\
public:							\
  absorb_##name() {flux_fun=&flux_fun_##name;};		\
};							\
							\
class bcconv_adv_##name : public BcconvAdv {		\
public:							\
  bcconv_adv_##name() {flux_fun=&flux_fun_##name;};	\
}

#define ADVECTIVE_ELEMSET_FM2(name)			\
FluxFunctionFM2 flux_fun_##name;			\
							\
class volume_##name : public AdvectiveFM2 {		\
public:							\
  volume_##name() {flux_fun=&flux_fun_##name;};		\
};							\
class bcconv_adv_##name : public BcconvAdvFM2 {		\
public:							\
  bcconv_adv_##name() {flux_fun=&flux_fun_##name;};	\
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Add here declarations for further advective elemsets. 
/// Euler equations for inviscid (Gas dynamics eqs.)

#define ADVDIF_ELEMSET(name)				\
class name##_ff_t : public AdvDifFF {			\
public:							\
  int operator()(ADVDIFFF_ARGS);			\
};							\
							\
class advdif_##name : public AdvDif {			\
public:							\
  advdif_##name() {adv_diff_ff = new name##_ff_t;};	\
};							\
							\
class bcconv_adv_##name : public BcconvAdv {		\
public:							\
  bcconv_adv_##name() {adv_diff_ff = new name##_ff_t;};	\
};

ADVDIF_ELEMSET(advecfm2);	// linear advective diffusive 
ADVDIF_ELEMSET(burgers);	// 1D scalar Burgers equation
ADVDIF_ELEMSET(swfm2t);	        // shallow water turbulent

class wall_swfm2t : NewElemset { 
public: 
  NewAssembleFunction new_assemble;
  ASK_FUNCTION;
};

// Global parameters
struct GlobParam {
  double alpha,Dt;
};

class newadvecfm2_ff_t : public NewAdvDifFF {
private:  
  int shock_capturing,na,nd,nc;
  FastMat2 D_jac_l, C_jac_l, tmp0;
  double tau_fac;
  FastMat2 u,u2,Uintri,AA,Ucpy;
  vector<double> djacv,cjacv;
  double *djacvp,*cjacvp;
  ElementIterator element;
  Property advective_jacobians_prop;
  const double *advjac;
  int ndim,ndof;
public:
  class UPerField;
  friend class UPerField;
  class UPerField : public AJac {
    newadvecfm2_ff_t &ff;
    FastMat2 tmp;
  public:
    UPerField(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_flux,comp_A_grad_U,comp_A_grad_N;
  };
  newadvecfm2_ff_t(NewAdvDif *elemset_) : NewAdvDifFF(elemset_) {};
  void start_chunk(int ret_options);
  void element_hook(ElementIterator &element);
  void compute_flux(COMPUTE_FLUX_ARGS);
};

class newadvdif_advecfm2 : public NewAdvDif {
public:
  newadvdif_advecfm2() : NewAdvDif(new newadvecfm2_ff_t(this)) {};
  // newadvdif_advecfm2() {adv_diff_ff = new newadvecfm2_ff_t(this);};
};

void log_transf(FastMat2 &true_lstate,const FastMat2 &lstate,
		const int nlog_vars,const int *log_vars);


#endif
