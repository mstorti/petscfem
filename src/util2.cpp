/*
  This file belongs to he PETSc - FEM package a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/
  
#include <stdio.h>
#include <cassert>
#include <map>
#include <vector>
#include <malloc.h>

#ifdef RH60
#include "libretto.h"
#else
#include <libretto/libretto.h>
#endif
#include <libretto/darray.h>

#include <sles.h>
#include <newmatio.h>
#undef HAVE_MEMMOVE // para que no chille al incluir petsccfonf.h

extern "C" {
#include "machine.h"
#include <matrix.h>
#include <matrix2.h>
}

#include "texthash.h"
#include "fstack.h"
#include "util2.h"

#undef __FUNC__
#define __FUNC__ "mat_function(const Matrix &A,Matrix &funA,scalarfun)"
int mat_function(const Matrix &A,Matrix &funA,scalarfun fun,void * user_data) {

  int n = A.Nrows();
  Matrix Lambda(n,2), V(n,n), Vim(n,n);
  non_symm_eigenvals(A,Lambda,V,Vim);
  DiagonalMatrix D(n);

  for (int k=1; k<=n; k++) {
    assert (Lambda(k,2)==0);
    D(k,k) = fun(Lambda(k,1),user_data);
  }

  funA = V * D * V.i();

}

#undef __FUNC__
#define __FUNC__ "mat_function(const Matrix &A,Matrix &funA,vectorfun)"
int mat_function(const Matrix &A,Matrix &funA,vectorfun fun,void * user_data) {

  int n = A.Nrows();
  Matrix Lambda(n,2), V(n,n), Vim(n,n);
  non_symm_eigenvals(A,Lambda,V,Vim);
  ColumnVector v(n),fv(n);
  DiagonalMatrix D(n);

  for (int k=1; k<=n; k++) {
    assert (Lambda(k,2)==0);
    v(k) = Lambda(k,1);
  }

  fv = fun(v,user_data);
  
  for (int k=1; k<=n; k++) {
    D(k,k) = fv(k);
  }

  funA = V * D * V.i();

}

double fabs_1(double a, void* p) {
  return fabs(a);
}

int matrix_abs(const Matrix &A,Matrix &absA) {
  mat_function(A,absA,&fabs_1,NULL);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Computes eigenvalues of a generic (may-be non-symmetric) matrix. 
    @author M. Storti
    @param A (input) n x n matrix. 
    @param lambda (input) n x 2 matrix containing real and imaginary
    parts of the eigenvalues. 
    @param V (input) n x n change of base matrix. 
*/

// Repeated here because `fem.h' includes <vector> and this collides
// with `meschach'
#define SHV(x) cout << #x ": " << x << endl

#undef __FUNC__
#define __FUNC__ "non_symm_eigenvals(const Matrix &,Matrix &,Matrix &)"
int non_symm_eigenvals(const Matrix &A,Matrix &lambda,Matrix &Vre,
		       Matrix &Vim) {

  MAT	*AA, *Q, *X_re, *X_im;
  VEC	*evals_re, *evals_im;
  double tol = 1e-7;

  int m,n;
  m=A.Nrows();
  n=A.Ncols();
  assert(m==n);
  AA = m_get(m,m);
  Q = m_get(m,m);
  X_re = m_get(m,m);
  X_im = m_get(m,m);
  evals_re = v_get(m);
  evals_im = v_get(m);

  double ev_shift = 0.;
  
  while (1) {
    double *p = A.Store();
    for (int k=0; k< m; k++) {
      for (int l=0; l< m; l++) {
	AA->me[k][l] = *(p++);
      }
      AA->me[k][k] += ev_shift;
    }


    /* compute Schur form: A = Q*T*Q^T */
    schur(AA,Q);
 
    /* extract eigenvalues */
    schur_evals(AA,evals_re,evals_im);

    double min_la, max_la, kmax, absla;
    for (int k=0; k<m; k++) {
      absla = sqrt(evals_re->ve[k]*evals_re->ve[k] 
		   +evals_im->ve[k]*evals_im->ve[k]);
      if (k==0 || absla <min_la)
	min_la = absla;
      if (k==0 || absla >min_la) {
	kmax = k;
	max_la = absla;
      }
    }

    if (min_la < tol*max_la) {
      ev_shift = 2*max_la;
      continue;
    }

    schur_vecs(AA,Q,X_re,X_im);

    for (int k=0; k<m; k++) {
      lambda(k+1,1) = evals_re->ve[k] - ev_shift;
      lambda(k+1,2) = evals_im->ve[k];
      for (int l=0; l<m; l++) {
	Vre(k+1,l+1) = X_re->me[k][l];
	Vim(k+1,l+1) = X_im->me[k][l];
      }
    }
    if (min_la >= tol*max_la) break;
  }

  M_FREE(AA);
  M_FREE(Q);
  M_FREE(X_re);
  M_FREE(X_im);
  V_FREE(evals_re);
  V_FREE(evals_im);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int vector_divide(res,a_mass,&dx)"
int vector_divide(Vec &res,Vec a_mass) {
  
  int lsize, ierr;
  double *ares,*aa_mass,*adx;
  ierr = VecGetLocalSize(res,&lsize); CHKERRQ(ierr); 

  // localize vectors
  ierr = VecGetArray(res,&ares);
  ierr = VecGetArray(a_mass,&aa_mass);
  for (int k=0; k<lsize; k++) {
    ares[k] /= aa_mass[k];
  }
  
  ierr = VecRestoreArray(res,&ares);
  ierr = VecRestoreArray(a_mass,&aa_mass);

}

double pw4(double x) {
  double x2 = x*x;
  return x2*x2;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void show_mallinfo (char *s=NULL)"
void show_mallinfo (char *s=NULL) {
  if (s!=NULL) printf("%s\n",s);
  printf("========== MEM USE ================================\n"
	 "Total memory (sbrk) [arena]: %d b\n"
	 "Alloc with malloc [uordblks]: %d chunks\n"
	 "Occup. by `free' (not in use) [fordblks]: %d chunks\n"
	 "===================================================\n",
	 mallinfo().arena,
	 mallinfo().uordblks,
	 mallinfo().fordblks);
}

#if 0 // Definicion completa
void show_mallinfo () {
  printf("========== MEM USE ================================\n"
	 "Total memory (sbrk) [arena]: %d b\n"
	 "Not in use [ordblks]: %d chunks\n"
	 "Alloc with mmap [hblkhd]: %d b\n"
	 "Alloc with malloc [uordblks]: %d chunks\n"
	 "Occup. by `free' (not in use) [fordblks]: %d chunks\n"
	 "Top most releaseable chunk [keepcost]: %d b(?)\n"
	 "===================================================\n",
	 mallinfo().arena,mallinfo().ordblks,
	 mallinfo().hblkhd,mallinfo().uordblks,
	 mallinfo().fordblks,mallinfo().keepcost);
}
#endif


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void read_double_array(vector<double> ,char *)" 
void read_double_array(vector<double> & ajac,char * advje) {
  char *buf= new char[strlen(advje)+1];
  strcpy(buf,advje);
  char *token;
  double val;
  for (int k=0; ; k++) {
    token = strtok(k==0 ? buf : NULL ," ");
    if (token==NULL) break;
    int nread = sscanf(token,"%lf",&val);
    if (nread!=1) break;
    ajac.push_back(val);
  }
  delete[] buf;
}

