//__INSERT_LICENSE__
//$Id: sparse.cpp,v 1.8 2001/09/21 16:54:32 mstorti Exp $

#include "sparse.h"

#if 0
void SeqMat::set_value(int row,int col,Scalar value,InsertMode
		       mode=ADD_VALUES) {
  ColIt col_it;
  RowIt row_it = find(row);
  Row r;
  if (row_it == end()) {
    r.clear();
    r.insert(col_p(col,value));
    insert(row_p(row,r));
  } else {
    col_it = row_it->find(col);
    if (col_it == row_it->end()) {
      col_it->insert(col_p(col,value));
    } else {
      col_it->second += value;
    }
  }
}

void SeqMat::create(Darray *da,const Dofmap *dofmap_,int debug_compute_prof=0) {
  MPI_Comm_size (PETSC_COMM_WORLD, &size);
  assert(size==1);
  nrows = dofmap->neq;
  ncols = dofmap->neq;
}
#endif


namespace Sparse {

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

  double Mat::get(int j,int k) const {
    assert(j<nrows);
    const RowCIt J = find(j);
    if (J == end()) {
      return 0.;
    } else {
      return J->second.get(k);
    }
  }

  int Vec::empty() const {
    VecCIt i,e;
    
    e = end();
    for (i=begin(); i!=e; i++) 
      if (i->second != 0.) return 0;
    return 1;
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

    e = end();
    for (i=begin(); i!=e; i++) {
      printf("row %d -> ",i->first);
      i->second.print_g(0,NULL,": ","  ","\n");
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

}

/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "/u/mstorti/PETSC/petscfem/tools/pfcpp")
  End: 
*/
