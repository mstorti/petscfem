//__INSERT_LICENSE__
//$Id: sparse.cpp,v 1.2 2001/09/20 20:13:24 mstorti Exp $

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
  Vec::Vec(const Indx &I,const Vec &v) {
    assert(0); // code here
  };

  double Vec::get(int j) const {
    const VecCIt J = find(j);
    if (J == end()) {
      return 0.;
    } else {
      return J->second;
    }
  }

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
      insert(VecP(j,v));
    } else {
      J->second = v;
    }
    return *this;
  }

  void Vec::print() {
    VecCIt J;
    printf("length: %d\n",len);
    for (J=begin(); J!=end(); J++) {
      printf("%d -> %f\n", J->first, J->second);
    }
  }

}

/*
  Local Variables: 
  eval: (setq c-macro-preprocessor "/u/mstorti/PETSC/petscfem-beta-1.93/tools/pfcpp")
  End: 
*/
