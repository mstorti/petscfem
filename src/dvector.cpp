//__INSERT_LICENSE__
//$Id: dvector.cpp,v 1.7 2005/01/17 23:51:25 mstorti Exp $
 
#include <src/dvector.h>
#include <src/dvector2.h>

// Explicit instantiation for `int' and `double'
template class dvector<double>;
template class dvector<int>;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dvector<double>::read(FILE *fid,double &t) {
  int nread = fscanf(fid,"%lf",&t);
  return nread!=1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dvector<int>::read(FILE *fid,int &t) {
  int nread = fscanf(fid,"%d",&t);
  return nread!=1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dvector<float>::read(FILE *fid,float &t) {
  int nread = fscanf(fid,"%f",&t);
  return nread!=1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dvector<double>::printe(FILE *fid,double t) {
  int ierr = fprintf(fid,"%.12g",t);
  return ierr<0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dvector<int>::printe(FILE *fid,int t) {
  int ierr = fprintf(fid,"%d",t);
  return ierr<0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dvector<float>::printe(FILE *fid,float t) {
  int ierr = fprintf(fid,"%.6g",t);
  return ierr<0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double dvector<double>::new_t_elem() { return 0.0; }

int dvector<int>::new_t_elem() {return 0; }

float dvector<float>::new_t_elem() { return 0.0; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
dvector<double>& dvector<double>::
axpy(double alpha,const dvector<double> &w) {
  assert(size()==w.size());
  int n = size();
  for (int j=0; j<n; j++) 
    ref(j) += alpha*w.ref(j);
  return *this;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
dvector<double>& dvector<double>::
scale(double alpha) {
  int n = size();
  for (int j=0; j<n; j++) 
    ref(j) *= alpha;
  return *this;
}

