// -*-mode: c++ -*-
//__INSERT_LICENSE__
//$Id: enthalpy.h,v 1.7 2005/02/23 01:40:34 mstorti Exp $
#ifndef PETSCFEM_ENTHALPY_H
#define PETSCFEM_ENTHALPY_H

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

/// User defined EF
class user_def_ef_t : public EnthalpyFun {
  /// Aux var. identity of size ndof
  FastMat2 eye_ndof,htmp1,htmp2;
  /// The actual Cp
  double Cp;
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
  // Return the Cp
  void get_Cp(FastMat2 &Cp_a);
};

#endif
