//__INSERT_LICENSE__
//$Id: iisdmat2.cpp,v 1.1 2003/08/31 21:21:09 mstorti Exp $
// fixme:= this may not work in all applications
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
extern int MY_RANK,SIZE;

#include "libretto.h"
#include <petscmat.h>

#include <src/fem.h>
#include <src/utils.h>
#include <src/dofmap.h>
#include <src/elemset.h>

#include <src/iisdmat.h>

int any_A_LL_other_stop = 0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "IISDMat::set_value_a"
int IISDMat::set_values_a(int nrows,int *idxr,int ncols,int *idxc,
			  PetscScalar *values, InsertMode mode) {
  int row_indx,col_indx,row_t,col_t;
  int nr[2], nc[2];

  // Should call a `isp_set_values()' here
  assert(nlay==0 || nlay==1);
  insert_mode = mode;  

  // Mapping of rows
  nr[L]=0;
  nr[I]=0;
  for (int jr=0; jr<nrows; jr++) {
    map_dof_fun(idxr[jr],row_t,row_indx);
    int jrl = nr[row_t]++;
    dvector<int> *indx = indxr[row_t];
    grow_mono<int>(*indx,nr[row_t]);
    indx->e(jrl) = row_indx;
    dvector<int> *jndx = jndxr[row_t];
    grow_mono<int>(*jndx,nr[row_t]);
    jndx->e(jrl) = jr;
  }

  // Mapping of columns
  nc[L]=0;
  nc[I]=0;
  for (int jc=0; jc<ncols; jc++) {
    map_dof_fun(idxc[jc],col_t,col_indx);
    int jcl = nc[col_t]++;
    dvector<int> *indx = indxc[col_t];
    grow_mono<int>(*indx,nc[col_t]);
    indx->e(jcl) = col_indx;
    dvector<int> *jndx = jndxc[col_t];
    grow_mono(*jndx,nc[col_t]);
    jndx->e(jcl) = jc;
  }

#if 0
  printf("nr[L] %d, nr[I] %d\n",nr[L],nr[I]);
  printf("indxr[L], jndxr[L]\n");
  for (int jr=0; jr<nr[L]; jr++) 
    printf("%d %d\n",indxr[L]->e(jr),jndxr[L]->e(jr));
  for (int jr=0; jr<nr[I]; jr++) 
    printf("%d %d\n",indxr[I]->e(jr),jndxr[I]->e(jr));

  printf("nc[L] %d, nc[I] %d\n",nc[L],nc[I]);
  printf("indxc[L], jndxc[L]\n");
  for (int jc=0; jc<nc[L]; jc++) 
    printf("%d %d\n",indxc[L]->e(jc),jndxc[L]->e(jc));
  for (int jc=0; jc<nc[I]; jc++) 
    printf("%d %d\n",indxc[I]->e(jc),jndxc[I]->e(jc));
  PetscFinalize();
  exit(0);
#endif

  // We consider first all blocks other than the LL one
  // This is considered aside
  for (row_t=0; row_t<2; row_t++) {
    for (col_t=0; col_t<2; col_t++) {
      if (col_t==L && row_t==L) continue;
      dvector<double> *vv = v[row_t][col_t];
      int nrr = nr[row_t];
      int ncc = nc[col_t];
      // Load values in matrix
      grow_mono(*vv,nr[row_t]*nc[col_t]);
      double *vvv = vv->buff();
      for (int jrl=0; jrl<nrr; jrl++) { 
	int jr = jndxr[row_t]->e(jrl);
	int *jndx = jndxc[col_t]->buff();
	double *values_jrow = values + jr*ncols;
	int *jc = jndx;
	int *jc_end = jc + ncc;
	while (jc<jc_end) *vvv++ = *(values_jrow + (*jc++));
      }
      // This is for debugging
      //#define LOAD_VALUES
#ifdef LOAD_VALUES
      ierr = MatSetValues(*(AA[row_t][col_t]),
			  nrr,indxr[row_t]->buff(),ncc,
			  indxc[col_t]->buff(),vv->buff(),mode);
      if (ierr) return ierr;
#endif
    }
  }    

  // Not implemented yet
  assert (local_solver != SuperLU);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // A_LL block

  // This will flag if we will consider the special case of having
  // LL elements to be sent to other processors via the A_LL_other
  // container
  int any_A_LL_other = 0;

  int nrr = nr[L];
  int ncc = nc[L];

  // Fix row matrix indices (block LL is shifted)
  int *p = indxr[L]->buff();
  int *pe = p + nrr;
  while (p<pe) {
    *p -= n_locp;
    if (*p < 0 || *p >= n_loc) any_A_LL_other = 1;
    p++;
  }

  // Fix column matrix indices (same reason)
  p = indxc[L]->buff();
  pe = p + ncc;
  while (p<pe) {
    *p -= n_locp;
    if (*p < 0 || *p >= n_loc) any_A_LL_other = 1;
    p++;
  }

  // Fix the matrices and send to A_LL_other if special case
  if (any_A_LL_other) {
    // send to A_LL_other
    any_A_LL_other_stop = 1;
    
    int *jrp = jndxr[L]->buff();
    int *indxrp = indxr[L]->buff();
    int *jcp = jndxc[L]->buff();
    int *indxcp = indxc[L]->buff();
    for (int jrl=0; jrl<nrr; jrl++) { 
      int jr = jrp[jrl];
      for (int jcl=0; jcl<ncc; jcl++) { 
	int jc = jcp[jcl];
	int indxr = indxrp[jrl];
	int indxc = indxcp[jcl];
	if (indxr <0 || indxr >=n_loc || indxc <0 || indxc >=n_loc) {
	  double w = values[jr*ncols+jc];
	  // A_LL_other->insert_val(indxr+n_locp,indxc+n_locp,w);
	  PetscSynchronizedPrintf(PETSC_COMM_WORLD,
				  "[%d] NEW, sending to A_LL_other %d,%d,%g\n",
				  MY_RANK,idxr[jr],idxc[jc],w);
	}
      }
    }
  }

  dvector<double> *vv = v[L][L];
  // Load values in matrix
  grow_mono(*vv,nr[L]*nc[L]);
  double *vvv = vv->buff();
  for (int jrl=0; jrl<nrr; jrl++) { 
    int jr = jndxr[L]->e(jrl);
    int *jndx = jndxc[L]->buff();
    double *values_jrow = values + jr*ncols;
    int *jc = jndx;
    int *jc_end = jc + ncc;
    while (jc<jc_end) *vvv++ = *(values_jrow + (*jc++));
  }

#ifdef LOAD_VALUES
  ierr = MatSetValues(*(AA[L][L]),
		      nrr,indxr[L]->buff(),ncc,
		      indxc[L]->buff(),vv->buff(),mode);
  if (ierr) return ierr;
#endif

  return 0;
}
