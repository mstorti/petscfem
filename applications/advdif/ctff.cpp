//__INSERT_LICENSE__
// $Id: ctff.cpp,v 1.1 2003/10/09 12:40:09 mstorti Exp $
/** Cutoff function. It is very near to ${\rm ctff(x)}\approx \rm tol$ for
    $x<0$ and ${\rm ctff}(x)=x$ for $x\gg \rm tol$. 
*/ 
double ctff(double x, double & diff_ctff, double tol) {
  double r=x/tol-1.; 
  double ee,vaux,ret;
  if (fabs(r)<1e-7) {
    ret = (1.+0.5*exp(r)/(1+(1./6.)*r*r))*tol;
    ee  = tol*tol;
//    vaux = 7.0*ee+x*x-2.0*x*tol;
//    diff_ctff = 3.0*ee*exp(r)*(9.0*ee+x*x-4.0*x*tol)/vaux/vaux;
    vaux = (1.+1./6.*r*r);
    diff_ctff  = 0.5*exp(r)*(vaux-1./6.*2*r)/vaux/vaux;
    // dfx1dx  = 0.5*exp(r).*(1+1/6*r.^2-1/6*2*r)./(1+1/6*r.^2).^2;
  } else if (r>0) {
    ret =  (x-tol)/(1.-exp(-2.*r))+tol;
    vaux  = exp(-2.*r); 
    diff_ctff = 1.0/(1.0-vaux)*(1.0-2.*r*vaux/(1.0-vaux));
    // dfx21dx = (1-2*r.*exp(-2*r)-exp(-2*r))./(1-exp(-2*r)).^2;
  } else {
    ee = exp(2.*r);
    ret = (x-tol)*ee/(ee-1.)+tol;
    vaux = ee-1.0;
    diff_ctff = ee/vaux*(1.0-2.0*r/vaux);
    // dfx22dx = (exp(2*r)./(exp(2*r)-1).^2).*(exp(2*r)-1-2*r);
    // printf("ctff(%g,%g) = %g\n",x,tol,ret);
  }
  return ret;
}
