// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: advective.h,v 1.3 2001/05/30 03:58:38 mstorti Exp $
 
#ifndef ADVECTIVE_H
#define ADVECTIVE_H


/** The jacobians of the flux functions. This is an array
    of ndim matrices of ndof x ndof entries each. 
*/
#define AJAC(jd) (*A_jac[(jd)-1])

/// Standard arguments for the typical flux function. Newmat version. 
#define FLUX_FUN_ARGS const RowVector &U,int ndim,const Matrix &iJaco, \
	      Matrix &H, Matrix &grad_H, \
	      Matrix &flux,vector<Matrix *> A_jac, \
	      Matrix &A_grad_U, Matrix &grad_U,  \
	      Matrix &G_source,Matrix &tau_supg, double &delta_sc, \
	      double &lam_max, \
	      TextHashTable *thash, \
              Matrix &nor,Matrix &lambda,Matrix &Vr,Matrix &Vr_inv, \
              double *propel, \
	      void *user_data,int options, int &start_chunk, int & ret_options

#define FLUX_FUN_ARGS_GENER(ML_) const ML_ &U,int ndim,const ML_ &iJaco, \
	      ML_ &H, ML_ &grad_H, \
	      ML_ &flux,vector<ML_ *> A_jac, \
	      ML_ &A_grad_U, ML_ &grad_U,  \
	      ML_ &G_source,ML_ &tau_supg, double &delta_sc, \
	      double &lam_max, \
	      TextHashTable *thash, \
              ML_ &nor,ML_ &lambda,ML_ &Vr,ML_ &Vr_inv, \
              double *propel, \
	      void *user_data,int options,int &start_chunk, int & ret_options

#define FLUX_FUN_ARGS_FM2 const  FastMat2  &U, \
	      int ndim, const  FastMat2  &iJaco,  FastMat2  &H, \
	      FastMat2  &grad_H,  FastMat2  &flux, \
	      FastMat2  &A_jac,  FastMat2  &A_grad_U, \
	      FastMat2  &grad_U,  FastMat2  &G_source, \
	      FastMat2  &tau_supg, double &delta_sc, double &lam_max, \
	      TextHashTable *thash,  FastMat2  &nor, FastMat2  &lambda, \
	      FastMat2  &Vr, FastMat2  &Vr_inv, double *propel, \
	      void *user_data,int options,int &start_chunk, int & ret_options

#define FLUX_FUN_ARGS_FM FLUX_FUN_ARGS_GENER(FastMat)

/// Standard arguments for the typical flux function (in the function call).
#define FLUX_FUN_CALL_ARGS_GENER(ML_) U##ML_,ndim, \
              iJaco##ML_,H##ML_,grad_H##ML_,flux##ML_,A_jac##ML_,A_grad_U##ML_, \
	      grad_U##ML_,G_source##ML_,tau_supg##ML_,delta_sc,lam_max,thash, \
              nor##ML_,lambda##ML_,Vr##ML_,Vr_inv##ML_,propel,user_data, \
	      options,start_chunk,ret_options

#define FLUX_FUN_CALL_ARGS_FM FLUX_FUN_CALL_ARGS_GENER(_FM)

#define FLUX_FUN_CALL_ARGS_FM2 FLUX_FUN_CALL_ARGS_GENER(_FM2)

#define FLUX_FUN_CALL_ARGS U,ndim,iJaco,H,grad_H,flux,A_jac,A_grad_U, \
	      grad_U,G_source,tau_supg,delta_sc,lam_max,thash, \
              nor,lambda,Vr,Vr_inv,propel,user_data, \
	      options,start_chunk,ret_options
	      
/// Type of flux function. 
typedef int FluxFunction (FLUX_FUN_ARGS);

typedef int FluxFunctionFM2 (FLUX_FUN_ARGS_FM2);

/** The class Advective is an Elemset class plus a
    pointer to a flux function. 
*/
class Advective : public Elemset { 
public:
  ASSEMBLE_FUNCTION;
  ASK_FUNCTION;
  FluxFunction *flux_fun;
};

class AdvectiveFM2 : public Elemset { 
public:
  ASSEMBLE_FUNCTION;
  ASK_FUNCTION;
  FluxFunctionFM2 *flux_fun;
};

class BcconvAdvFM2 : public Elemset { 
public:
  ASSEMBLE_FUNCTION;
  ASK_FUNCTION;
  FluxFunctionFM2 *flux_fun;
};

/** The class Absorb is an Elemset class plus a
    pointer to a flux function. 
*/
class Absorb : public Elemset { 
public:
  ASSEMBLE_FUNCTION;
  ASK_FUNCTION;
  FluxFunction *flux_fun;
};

class BcconvAdv : public Elemset { 
public:
  ASSEMBLE_FUNCTION;
  ASK_FUNCTION;
  FluxFunction *flux_fun;
};

enum flux_fun_opt {
  DEFAULT     = 0x0000,
  COMP_SOURCE = 0x0001,
  COMP_UPWIND = 0x0002,
  COMP_EIGENV = 0x0004,
  SCALAR_TAU  = 0x0008,
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define ADVECTIVE_ELEMSET(name) \
FluxFunction flux_fun_##name; \
 \
class volume_##name : public Advective { \
public: \
  volume_##name() {flux_fun=&flux_fun_##name;}; \
}; \
 \
class absorb_##name : public Absorb { \
public: \
  absorb_##name() {flux_fun=&flux_fun_##name;}; \
}; \
 \
class bcconv_adv_##name : public BcconvAdv { \
public: \
  bcconv_adv_##name() {flux_fun=&flux_fun_##name;}; \
}

#define ADVECTIVE_ELEMSET_FM2(name) \
FluxFunctionFM2 flux_fun_##name; \
 \
class volume_##name : public AdvectiveFM2 { \
public: \
  volume_##name() {flux_fun=&flux_fun_##name;}; \
}; \
class bcconv_adv_##name : public BcconvAdvFM2 { \
public: \
  bcconv_adv_##name() {flux_fun=&flux_fun_##name;}; \
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Add here declarations for further advective elemsets. 
/// Euler equations for inviscid (Gas dynamics eqs.)
ADVECTIVE_ELEMSET(euler);    
ADVECTIVE_ELEMSET_FM2(eulerfm2);    // con FastMat2
/// Shallow water equations. 
ADVECTIVE_ELEMSET(shallow);
ADVECTIVE_ELEMSET_FM2(shallowfm2);    // con FastMat2
ADVECTIVE_ELEMSET_FM2(shallowfm2t);    // con FastMat2 / turbulento
/// Advection of multiple scalar fields with a velocity field. 
ADVECTIVE_ELEMSET(advec);

// use FastMat2 cache or not in advecfm2/eulerfm2.cpp
#define USE_FASTMAT2_CACHE

#endif

