//__INSERT_LICENSE__
//$Id: fm2funm.h,v 1.4 2007/01/30 19:03:44 mstorti Exp $
#ifndef PETSCFEM_FM2FUNM_H
#define PETSCFEM_FM2FUNM_H

#define FastMat2_fund    FastMat2_fund_advdif
#define FastMat2_funm    FastMat2_funm_advdif
#define FastMat2_funm_ff FastMat2_funm_ff_advdif

class FastMat2_fund {
public:
  virtual ~FastMat2_fund() {}
private:
  int m;
  FastMat2 V,D,Vi,tmp,lambda,flambda;
public:
  void init(FastMat2 &A);
  void apply(const FastMat2 &A,FastMat2 &fA);
  void apply_ns(const FastMat2 &A,FastMat2 &fA);
  virtual void f(const FastMat2 &D,FastMat2 &fD)=0;
};

double FastMat2_funm_ff(double x,void *a);

class FastMat2_funm : public FastMat2_fund {
private:
  void f(const FastMat2 &D,FastMat2 &fD);
public:
  virtual double f(double)=0;
};


#endif
