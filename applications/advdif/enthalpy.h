// -*-mode: c++ -*-
#ifndef ENTHALPY_H
#define ENTHALPY_H

// Constant Cp, the same for all fields
class ScalarPerFieldEF : public EnthalpyFun {
  FastMat2 Cp,htmp1,htmp2;
  int ndof;
public:
  void init(int ndim,int ndof,int nel);
  void update(const double *ejac);
  void enthalpy(FastMat2 &H, FastMat2 &U);
  void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
		   double w);
  void comp_P_Cp(FastMat2 &P_Cp,FastMat2 &P_supg);
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
  void comp_P_Cp(FastMat2 &P_Cp,FastMat2 &P_supg);
};

#endif
