//__INSERT_LICENSE__
//$Id: sparse.cpp,v 1.15 2001/09/23 15:58:17 mstorti Exp $

#include "sparse.h"

using namespace Random;

namespace Sparse {

  GenVec::~GenVec() {};

  double GenVec::not_represented_val = 0.;
  double Mat::not_represented_val = 0.;
  
  ScalarFunObj::~ScalarFunObj() {};

  Scale scale_fun_obj;

  ScalarFunWrapper scalar_fun_wrapper;

  BinAssoc::~BinAssoc() {};

  Sum sum_bin_assoc;
  SumAbs sum_abs_bin_assoc;
//    SumSq sum_sq_bin_assoc;
//    SumPow sum_pow_bin_assoc;
  Max max_bin_assoc;
  MaxAbs max_abs_bin_assoc;
  Min min_bin_assoc;

  Accumulator::~Accumulator() {};
  SumSq sum_sq_accum;
  SumPow sum_pow_accum;

  double ScalarFunWrapper::fun(double v) const {
    if (sfd) {
      return (*sfd)(v,user_data);
    } else {
      return (*sf)(v);
    }
  } 

  void GenVec::print_f(const char *s = NULL) {
    if (s) printf("%s\n",s);
    int m = length();
    for (int j=0; j<m; j++) {
      printf("%d %f\n",j,get(j));
    }
  }

  /// print elements
  void GenVec::print(const char *s = NULL) {
    double v;
    if (s) printf("%s\n",s);
    int m = length();
    for (int j=0; j<m; j++) {
      v = get(j);
      if (v!=0.) printf("%d %f\n",j,v);
    }
  }
    
  GenVec & GenVec::prod(const Mat & a,const GenVec & v) {
    int m,n;
    RowCIt i,e;
    Vec w;

    m = a.rows();
    n = a.cols();
    assert(length()==m);
    assert(v.length()==n);

    e = a.end();
    for (i=a.begin(); i!=e; i++) {
      w = i->second;
      w.resize(n);
      set(i->first,w.dot(v));
    }

    return *this;
  }

  Indx::Indx(int m,int n,int k=1) {
    int j,v;
    assert(m>=0);
    assert(n>=0);
    // assert(k*(n-m)>0);
    if (n>=m) resize((n-m)/k+1);
    v = m;
    j = 0;
    while (k*(v-n) <= 0) {
      (*this)[j++] = v;
      v += k;
    }
  }

#if 0
  Vec::Vec(const Indx &I,const Vec &v) {
    assert(0); // code here
  };
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Vec::get"
  double Vec::get(int j) const {
    assert(j<len);
    return get_nc(j);
  }

  double Vec::get_nc(int j) const {
    const VecCIt J = find(j);
    if (J == end()) {
      return 0.;
    } else {
      return J->second;
    }
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "Vec::print"
  void Vec::print(const char * s = NULL) {
    print_g(0,s," -> ","  ","\n");
  }
#endif

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Vec::print_g"
  void Vec::print_g(int l,const char * s,const char * psep, const char * isep,
		  const char * lsep) const {
    VecCIt J;
    if (s) printf("%s\n",s);
    if (l) printf("length: %d\n",len);
    for (J=begin(); J!=end(); J++) {
      printf("%d%s%f%s", J->first, psep, J->second,isep);
    }
    printf("%s",lsep);
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Vec::print"
  void Vec::print(const char * s = NULL) {
    print_g(1,s," -> ","\n","-- end vec --\n");
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Vec::resize"
  Vec & Vec::resize(int n) {
    VecIt i;
    assert(n>=0);
    for (i=begin(); i!= end(); i++) 
      if (i->first >= n) erase(i);
    len = n;
    return *this;
  }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Vec::set"
  Vec & Vec::set(int j,double v) {
    VecIt J = find(j);
    if (j >= len ) {
      if (!grow_m) {
	printf("sparse: index %d out of range, length %d\n",j,len);
	abort();
      }
      len = j+1;
    }
    if (J == end()) {
      if (v!=0) insert(VecP(j,v));
    } else {
      if (v==0) {
	erase(J);
      } else {
	J->second = v;
      }
    }
    return *this;
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Vec::set"
  Vec & Vec::set(const Indx &I,const Vec & v) {
    int j;
    for (j=0; j<I.size(); j++) {
      set(I[j],v.get(j));
    }
    return *this;
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Vec::set"
  Vec & Vec::set(const Vec &v,const Indx &I)  {
    int j,m;
    clear();
    m = I.size();
    resize(m);
    for (j=0; j<m; j++) {
      set(j,v.get(I[j]));
    }
    return *this;
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ ""
  Vec & Vec::set(const Indx & I,const Vec & v,const Indx & K) {
    Vec tmp;
    tmp.set(v,K);
    set(I,tmp);
    return *this;
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ ""
  Vec & Vec::scale(double c) {
    VecIt i,e;
    if (c==0.) {
      clear();
    } else {
      e = end();
      for (i=begin(); i!=e; i++) 
	i->second *= c;
    }
    return *this;
  }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Vec::axpy"
  Vec & Vec::axpy(double c,const Indx & I,double a,const Vec & v) {
    int j,m,k;
    m = I.size();
    for (j=0; j<m; j++) {
      k = I[j];
      set(k,c*get(k)+a*v.get(j));
    }
    return *this;
  }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "Vec::set"
  Vec & Vec::axpy(double a,const Vec & v) {
    VecCIt j,e;
    double r;
    int k;
    assert(len == v.len);
    e = v.end();
    for (j=v.begin(); j!=e; j++) {
      k = j->first;
      r = get(k);
      set(k, r + a * j->second);
    }
    return *this;
  }

  int Vec::empty() const {
    VecCIt i,e;
    
    e = end();
    for (i=begin(); i!=e; i++) 
      if (i->second != 0.) return 0;
    return 1;
  }

  Vec & Vec::setc(const Mat &a,int j) {
    RowCIt r,e;

    assert (j < a.ncols ) ;

    clear();
    resize(a.ncols);
    // Clear the column
    e = a.end();
    for (r = a.begin(); r!=e; r++) 
      set(r->first,r->second.get(j));
    return *this;
  }

  Vec & Vec::set(const Mat & a,int j) {
    RowCIt it;
    clear();
    it = a.find(j);
    if (it!=a.end()) 
      copy(it->second).resize(a.ncols);
    return *this;
  }

#if 1
  double TOL;
  int purge_fun(pair<const int,double> &p) {
    // double v = fabs(p.second);
    // printf("v %f, TOL %f, v<tol %d\n",v,TOL,v<TOL);
    return fabs(p.second)<TOL;
  }
#endif

  Vec & Vec::purge(double tol = 1e-10) {
    VecIt i,e;
#if 0 // This doesn't compile. I think there is some problem
      // with the remove_if algorithm with maps
    i = remove_if(begin(),end(),&purge_fun);
    erase(i,end());
#elsif 0
    e = end();
    for (i=begin(); i!=end(); i++) {
      printf("second %f, erase ? %d\n",i->second,fabs(i->second)<tol);
      if (fabs(i->second) < tol) erase(i);
    }
#else
    VecIt j;
    TOL = tol;
    // Advance to the first remaining item
    while (1) {
      i = begin();
      if (!purge_fun(*i)) break;
      erase(i);
    }
    if (i!=end()) {
      while (1) {
	// look ahead at following item
	j = i; j++;
	if (j==end()) break;
	if (purge_fun(*j)) {
	  // erase it ...
	  erase(j);
	} else {
	  // ... or advance iterator 
	  i=j;
	}
      }
    }
#endif
    return *this;
  }

  Vec & Vec::random_fill(double fill=0.1,Generator & g=uniform) {
    int j,k;

    clear();
    for (j=0; j<fill*len; j++) {
      k = irand(0,len-1);
      set(k,g.get());
    }
    return *this;
  }

  double Vec::dot(const GenVec & w) const {
    VecCIt i,e;
    double sum;

    assert(len == w.length());

    // For efficiency we loop over the elements
    // of the vector that has less elements (w1)
    sum = 0.;
    e = end();
    for (i = begin(); i!=e; i++) 
      sum += (i->second) * w.get(i->first);
    return sum;
  }
  
  double Vec::dot(const Vec & w) const {
    VecCIt i,e;
    const Vec *w1,*w2;
    double sum;

    assert(len == w.len);

    // For efficiency we loop over the elements
    // of the vector that has less elements (w1)
    sum = 0.;
    if (size()>w.size()) {
      w1 = &w;
      w2 = this;
    } else {
      w1 = this;
      w2 = &w;
    }
    e = w1->end();
    for (i = w1->begin(); i!=e; i++) 
      sum += (i->second) * w2->get(i->first);
    return sum;
  }
  
  double Mat::get(int j,int k) const {
    assert(j<nrows);
    assert(k<ncols);
    const RowCIt J = find(j);
    if (J == end()) {
      return 0.;
    } else {
      return J->second.get_nc(k);
    }
  }

  /// Set element at position j
  Mat & Mat::set(int j,int k,double v) {
    RowIt J;
    Vec row;
    if (j >= nrows ) {
      if (!grow_m) {
	printf("sparse: index %d out of range, rows  %d\n",j,nrows);
	abort();
      }
      nrows = j+1;
    }

    if (k >= ncols ) {
      if (!grow_m) {
	printf("sparse: index %d out of range, cols  %d\n",j,nrows);
	abort();
      }
      ncols = k+1;
    }

    J = find(j);
    if (J == end()) {
      if (v!=0) {
	row.set(k,v);
	insert(RowP(j,row));
      }
    } else {
      J->second.set(k,v);
    }
    return *this;
  }

  void Mat::print(const char *s = NULL) {
    RowCIt i,e;
    if (s) printf("%s\n",s);

    printf("size %d %d \n",nrows,ncols);
    e = end();
    for (i=begin(); i!=e; i++) {
      printf("row %d -> ",i->first);
      i->second.print_g(0,NULL,": ","  ","\n");
    }
    printf("-- end mat --\n");
  }

  void Mat::print_f(const char *s = NULL) {
    RowCIt i,e;
    int j,k;
    if (s) printf("%s\n",s);

    printf("size %d %d \n",nrows,ncols);
    for (j=0; j<nrows; j++) {
      printf("row %d: ",j);
      for (k=0; k<ncols; k++) 
	printf("%f ",get(j,k));
      printf("\n");
    }
    printf("-- end mat --\n");
  }

  /// Resize vectors, truncates elements if greater than this value
  Mat & Mat::resize(int m,int n) {
    RowIt i,e;

    i = lower_bound(m);
    erase(i,end());
    
    e = end();
    for (i=begin(); i!=e; i++) i->second.resize(n);
    nrows = m;
    ncols = n;
    return *this;
  }

  Mat & Mat::setr(int j,Vec &v) {
    RowIt J = find(j);
    int v_empty;

    if (j >= nrows ) {
      if (!grow_m) {
	printf("sparse: index %d out of range, length %d\n",j,nrows);
	abort();
      }
      nrows = j+1;
    }
    v_empty = v.empty();

    if (J == end()) {
      if (!v_empty) insert(RowP(j,v));
    } else {
      if (v_empty) {
	erase(J);
      } else {
	J->second = v;
      }
    }
    return *this;
  }

  Mat & Mat::setc(int j,const Vec &v) {
    RowIt r;
    VecCIt k,e;

    if (j >= ncols ) {
      if (!grow_m) {
	printf("sparse: index %d out of range, ncols %d\n",j,ncols);
	abort();
      }
      ncols = j+1;
    }

    // Clear the column
    for (r=begin(); r!=end(); r++) {
      r->second.set(j,0.);
      if (r->second.empty()) erase(r);
    }

    // Set the column
    e = v.end();
    for (k=v.begin(); k!=e; k++) 
      set(k->first,j,k->second);

    return *this;
  }

  void Mat::getr(int j,Vec & v) const {
    RowCIt it;
    v.clear();
    it = find(j);
    if (it!=end()) 
      v.copy(it->second).resize(ncols);
  }

  Mat & Mat::setr(const Mat & a,Indx &J) {
    int j,m,p;
    Vec row;
    
    m = J.size();
    clear();
    resize(m,a.ncols);
    for (j=0; j<m; j++) {
      p = J[j];
      assert(p<a.nrows);
      a.getr(p,row);
      setr(j,row);
    }
    return *this;
  }

  Mat & Mat::setc(const Mat & a,Indx &J) {
    RowCIt r,e;
    int m;
    Vec row;

    m=J.size();
    clear();
    resize(a.nrows,m);

    e = a.end();
    for (r=a.begin(); r!=e; r++) {
      row.set(r->second,J);
      insert(RowP(r->first,row));
    }
    return *this;
  }

  Mat & Mat::id(double a=1.) {
    assert(nrows==ncols);
    clear();
    for (int j=0; j<nrows; j++) set(j,j,a);
    return *this;
  }

  Mat & Mat::diag(const Vec & v) {
    int m,j;
    VecCIt i,e;

    clear();
    m = v.length();
    resize(m,m);
    e = v.end();
    for (i=v.begin(); i!=e; i++) {
      j = i->first;
      set(j,j,i->second);
    }
    return *this;
  }

  int Mat::empty() const {
    RowCIt i,e;
    
    e = end();
    for (i=begin(); i!=e; i++) 
      if (!(i->second.empty())) return 0;
    return 1;
  }

  int Mat::size() const {
    RowCIt r,e;
    int s=0;
    
    e = end();
    for (r=begin(); r!=end(); r++) 
      s += r->second.size();
    return s;
  }

  Mat & Mat::random_fill(double fill=0.1,Generator & g=uniform) {
    int j,k,l,m,n;

    clear();
    m = rows();
    n = cols();
    for (j=0; j<fill*m*n; j++) {
      k = irand(0,m-1);
      l = irand(0,n-1);
      set(k,l,g.get());
    }
    return *this;
  }

  //Computes w += c * a * b
  Mat & Mat::prod(const Mat & a, const Mat & b,double c=1.) {
    RowCIt i,e;
    VecCIt jk,ee;
    Vec a_row,b_row,row;
    int j,k,n,m,p;
    double a_jk;

    m = rows();
    n = cols();
    p = a.cols();
    assert(m==a.rows());
    assert(n==b.cols());
    assert(p==b.rows());

    e = a.end();
    for (i=a.begin(); i!=e; i++) {
      j = i->first;
      a_row = i->second;
      ee = a_row.end();
      for (jk=a_row.begin(); jk!=ee; jk++) {
	k = jk->first;
	a_jk = jk->second;
	b.getr(k,b_row);
	b_row.resize(n);
	getr(j,row);
	row.resize(n);
	row.axpy(c*a_jk,b_row);
	setr(j,row);
      }
    }
    return *this;
  }

  Vec & Vec::apply(const ScalarFunObj & fun) {
    VecIt i;
    for (i=begin(); i!=end(); i++) 
      i->second = fun.fun(i->second);
    return *this;
  }

  /// Apply a scalar function with args to all elements
  Vec & Vec::apply(ScalarFunD *fun,void * user_data = NULL) {
    scalar_fun_wrapper.user_data = user_data;
    scalar_fun_wrapper.sf = NULL;
    scalar_fun_wrapper.sfd = fun;
    return apply(scalar_fun_wrapper);
  }

  /// Apply a scalar function with no args to all elements
  Vec & Vec::apply(ScalarFun *fun) {
    scalar_fun_wrapper.user_data = NULL;
    scalar_fun_wrapper.sf = fun;
    scalar_fun_wrapper.sfd = NULL;
    return apply(scalar_fun_wrapper);
  }

  double Vec::assoc(BinAssoc & op) const {
    double val;
    VecCIt i,e;
    int j,m;

    val = op.null;
    e = end();
    for (i=begin(); i!=e; i++) 
      val = op.op(val,i->second);

    // If some value is not in the map, then it represents
    // `not_represented_val' (probably 0).
    m = length();
    for (j=0; j<length(); j++) {
      if (find(j)!=e) {
	val = op.op(val,not_represented_val);
	break;
      }
    }

    return val;
  }

  void Vec::accum(double &v,Accumulator & acc) const {
    VecCIt i,e;
    int j,m;

    e = end();
    for (i=begin(); i!=e; i++) 
      acc.accum(v,i->second);

    // If some value is not in the map, then it represents
    // `not_represented_val' (probably 0).
    m = length();
    for (j=0; j<length(); j++) {
      if (find(j)!=e) {
	acc.accum(v,not_represented_val);
	break;
      }
    }

  }

  Mat & Mat::apply(const ScalarFunObj & fun) {
    RowIt i;
    for (i=begin(); i!=end(); i++) 
      i->second.apply(fun);
    return *this;
  }

  /// Apply a scalar function with args to all elements
  Mat & Mat::apply(ScalarFunD *fun,void * user_data = NULL) {
    scalar_fun_wrapper.user_data = user_data;
    scalar_fun_wrapper.sf = NULL;
    scalar_fun_wrapper.sfd = fun;
    return apply(scalar_fun_wrapper);
  }

  /// Apply a scalar function with no args to all elements
  Mat & Mat::apply(ScalarFun *fun) {
    scalar_fun_wrapper.user_data = NULL;
    scalar_fun_wrapper.sf = fun;
    scalar_fun_wrapper.sfd = NULL;
    return apply(scalar_fun_wrapper);
  }

  double Mat::assoc(BinAssoc & op) const {
    double val;
    RowCIt i,e;
    Vec row;
    int j,m,n;
    double w;

    n = cols();
    e = end();
    val = op.null;
    for (i=begin(); i!=e; i++) {
      row = i->second;
      row.resize(n);
      w = row.assoc(op);
      val = op.op(val,w);
    }

    // If some value is not in the map, then it represents
    // `not_represented_val' (probably 0).
    m = rows();
    for (j=0; j<m; j++) {
      if (find(j)!=e) {
	val = op.op(val,not_represented_val);
	break;
      }
    }

    return val;
  }

  void Mat::accum(double &v,Accumulator & acc) const {
    RowCIt i,e;
    Vec row;
    int j,m,n;
    double w;

    n = cols();
    e = end();
    for (i=begin(); i!=e; i++) {
      row = i->second;
      row.resize(n);
      row.accum(v,acc);
    }

    // If some value is not in the map, then it represents
    // `not_represented_val' (probably 0).
    m = rows();
    for (j=0; j<m; j++) {
      if (find(j)!=e) {
	acc.accum(v,not_represented_val);
	break;
      }
    }

  }
}

/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "/u/mstorti/PETSC/petscfem/tools/pfcpp")
  End: 
*/
