// -*-mode: c++ -*-
//__INSERT_LICENSE__
//$Id: enthalpy.h,v 1.6 2003/01/08 15:54:25 mstorti Exp $
#ifndef ENTHALPY_H
#define ENTHALPY_H

/// Linear Constant Cp, the same for all fields
class ScalarPerFieldEF : public EnthalpyFun {
  /// Stores a cp for each field
  FastMat2 Cp,htmp1,htmp2;
  /// the number of fields
  int ndof;
public:
  /// initializes the object
  void init(int ndim,int ndof,int nel);
  /// sets the Cp coefficients for each field
  void update(const double *ejac);
  /// Scales each $U$ value by the corresponding $C_p$ coefficient.
  void enthalpy(FastMat2 &H);
  /// Efficient implementation
  void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w);
  /// Efficient implementation
  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);
};

/// A general Cp matrix. See base class documentation.
class FullEF : public EnthalpyFun {
  FastMat2 Cp,htmp1,htmp2;
public:
  /// Initializes
  void init(int ndim,int ndof,int nel);
  /// Full implementation
  void update(const double *ejac);
  /// Full implementation
  void enthalpy(FastMat2 &H);
  /// Full implementation
  void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w);
  /// Full implementation
  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);
};

#endif
