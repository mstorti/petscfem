// -*-mode: c++ -*-
#ifndef ENTHALPY_H
#define ENTHALPY_H

class EnthalpyFun {
public:
  virtual void update(const double *) {};
  virtual void enthalpy(FastMat2 &H, FastMat2 &U)=0;
  virtual void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
			   double w)=0;
};

// Constant Cp for all the fields
class GlobalScalarEF : public EnthalpyFun {
  FastMat2 eye_ndof,htmp1,htmp2;
  double Cp;
public:
  void init(int ndim,int ndof,int nel,double Cp=1.);
  void update(const double *Cp_) {Cp=*Cp_;};
  void enthalpy(FastMat2 &H, FastMat2 &U);
  void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
			   double w);
};

// Constant Cp=1 for all the fields. Identity relation between H and T
class IdentityEF : public GlobalScalarEF {
public:
  void update(const double *Cp_) {};
  void enthalpy(FastMat2 &H, FastMat2 &U) {H.set(U);};
};

// Constant Cp, the same for all fields
class ScalarPerFieldEF : public EnthalpyFun {
  FastMat2 Cp,htmp1,htmp2;
public:
  void init(int ndim,int ndof,int nel);
  void update(const double *ejac);
  void enthalpy(FastMat2 &H, FastMat2 &U);
  void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
		   double w);
};

// A general Cp matrix 
class FullEF : public EnthalpyFun {
  FastMat2 Cp,htmp1,htmp2;
public:
  void init(int ndim,int ndof,int nel);
  void update(const double *ejac);
  void enthalpy(FastMat2 &H, FastMat2 &U);
  void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
		   double w);
};

#endif
