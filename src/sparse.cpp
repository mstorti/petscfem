//__INSERT_LICENSE__
//$Id: sparse.cpp,v 1.6 2001/09/21 02:23:35 mstorti Exp $

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

  Vec::Vec(const Indx &I,const Vec &v) {
    assert(0); // code here
  };

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
#undef __FUNC__
#define __FUNC__ "Vec::print"
  void Vec::print(char * s = NULL) {
    VecCIt J;
    if (s) printf("%s\n",s);
    printf("length: %d\n",len);
    for (J=begin(); J!=end(); J++) {
      printf("%d -> %f\n", J->first, J->second);
    }
    printf("--\n");
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
#define __FUNC__ "Vec::set"
  Vec & Vec::axpy(double c,const Indx & I,double a,const Vec & v) {
    int j,m,k;
    m = I.size();
    for (j=0; j<m; j++) {
      k = I[j];
      set(k,c*get(k)+a*v.get(j));
    }
    return *this;
  }

  Vec & Vec::axpy(double a,const Vec & v) {
    VecCIt j,e;
    int k;
    assert(len == v.len);
    e = v.end();
    for (j=v.begin(); j!=e; j++) {
      k = j->first;
      set(k,get(k)+a*j->second);
    }
    return *this;
  }

#if 0
  Vec & Vec::axpy(int j,double a,double v) {
    double o;
    set(j,get(j)+a*
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
#endif
}

/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "/u/mstorti/PETSC/petscfem-beta-1.93/tools/pfcpp")
  End: 
*/
