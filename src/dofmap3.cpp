//__INSERT_LICENSE__
//$Id: dofmap3.cpp,v 1.5 2003/08/31 13:23:56 mstorti Exp $

#include <cassert>

using namespace std;

#include "fem.h"
#include "getprop.h"
#include "fstack.h"
#include "dofmap.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Dofmap::freeze() {

  // Initialization
  idmap2_dv.mono(nnod*ndof);
  idmap2 = idmap2_dv.buff();
  special_ptr_dv.clear();
  sp_eq_dv.clear();
  coefs_dv.clear();
  one_coef = 1.0;

  IdMapRow row;
  // Posicion in `sp_eq' and `coefs'
  int last = 0;
  // Index of special node/field
  int sp_indx=0;
  for (int node=1; node<=nnod; node++) {
    for (int kdof=1; kdof<=ndof; kdof++) {
      get_row(node,kdof,row);
      int s = row.size();
#if 0
      printf("%d/%d: ",node,kdof);
      for (int j=0; j<s; j++) {
	IdMapEntry &e = row[j];
	assert(e.j>=1 && e.j<=neqtot);
	printf("%d %f, ",e.j,e.coef);
      }
      printf("\n");
#endif
      int indx = edof(node,kdof)-1;
      if (s==1) {
	IdMapEntry &e = row[0];
	if (e.coef==1.0) {
	  // node/field combination is regular (mapped
	  // with coef==1.0 to a single equation)
	  idmap2[indx] = e.j;
	  continue;
	}
      }
      // node/field is special
      idmap2[indx] = neqtot + 1 + sp_indx++;
      special_ptr_dv.push(last);
      for (int j=0; j<s; j++) {
	IdMapEntry &e = row[j];
	sp_eq_dv.push(e.j);
	coefs_dv.push(e.coef);
      }
      last += s;
    }
  }
  special_ptr_dv.push(last);
  
  special_ptr_dv.defrag();
  special_ptr = special_ptr_dv.buff();
  sp_eq_dv.defrag();
  sp_eq = sp_eq_dv.buff();
  coefs_dv.defrag();
  coefs = coefs_dv.buff();
#if 0
  for (int node=1; node<=nnod; node++) {
    for (int kdof=1; kdof<=ndof; kdof++) {
      int ndof, s;
      const int *dofs;
      const double *coef;
      get_row(node,kdof,s,&dofs,&coef);
      printf("%d/%d: ",node,kdof);
      for (int j=0; j<s; j++) {
	printf("%d %f, ",dofs[j],coef[j]);
      }
      printf("\n");
    }
  }
#endif
}
	
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 	       
void Dofmap::get_row(int node,int kdof,int &neqs,
		     const int **dofpp, const double **coefpp) const {
  // Old - slow
  //  int *jndx = idmap2 + edof(node,kdof)-1;
  // New fast
  int *dofp = idmap2 + (node-1)*ndof + kdof-1;
  int indx = *dofp;
  if (indx<=neqtot) {
    neqs = 1;
    *dofpp = dofp;
    *coefpp = &one_coef;
  } else {
    indx -= (neqtot+1);
    int n1 = special_ptr[indx];
    neqs = special_ptr[indx+1] - n1;
    *dofpp = sp_eq+n1;
    *coefpp = coefs+n1;
  }
}
