//__INSERT_LICENSE__
//$Id: fm2funm.h,v 1.3 2007/01/30 19:03:44 mstorti Exp $
#ifndef FM2FUNM_H
#define FM2FUNM_H

#define FastMat2_fund    FastMat2_fund_ns
#define FastMat2_funm    FastMat2_funm_ns
#define FastMat2_funm_ff FastMat2_funm_ff_ns

class FastMat2_fund {
 public:
  virtual ~FastMat2_fund() {}
private:
  int m;
  FastMat2 V,D,tmp,lambda,flambda;
public:
  void init(FastMat2 &A);
  void apply(const FastMat2 &A,FastMat2 &fA);
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
