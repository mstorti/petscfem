//__INSERT_LICENSE__
//$Id: dvector.cpp,v 1.4 2003/08/28 22:43:16 mstorti Exp $
 
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
int dvector<double>::print(FILE *fid,double t) {
  int ierr = fprintf(fid,"%.12g\n",t);
  return ierr<0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dvector<int>::print(FILE *fid,int t) {
  int ierr = fprintf(fid,"%d\n",t);
  return ierr>=0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int dvector<float>::print(FILE *fid,float t) {
  int ierr = fprintf(fid,"%.6g\n",t);
  return ierr>=0;
}

