//__INSERT_LICENSE__
//$Id: dofmap3.cpp,v 1.13.10.1 2007/02/19 20:23:56 mstorti Exp $

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
void Dofmap::
get_row(int node,int kdof,int &neqs,
	const int **dofpp, 
	const double **coefpp) const {
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Dofmap::
qxpy(double *y,double *x, double alpha) {
  int nrow = nnod*ndof;
  int ncol = neqtot;
  int m;
  for (int j=0; j<nrow; j++) {
    int kdof = j % ndof + 1;
    int node = j / ndof + 1;
    const int *dofs;
    const double *coef;
    get_row(node,kdof,m,&dofs,&coef);
    double v = 0.0;
    for (int l=0; l<m; l++) 
      v += coef[l]*x[dofs[l]-1];
    y[j] += alpha*v;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Dofmap::
qtxpy(double *y,double *x, double alpha) {
  int nrow = nnod*ndof;
  int ncol = neqtot;
  int m;
  for (int j=0; j<nrow; j++) {
    int kdof = j % ndof + 1;
    int node = (j / ndof) + 1;
    const int *dofs;
    const double *coef;
    get_row(node,kdof,m,&dofs,&coef);
    for (int l=0; l<m; l++)
      y[dofs[l]-1] += alpha*coef[l]*x[j];
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
static void 
printv(double *v,int n,const char *s=NULL) {
  if (s) printf("%s\n",s);
  for (int j=0; j<n; j++)
    printf("%f\n",v[j]);
}
static double
diff(double *v, double *w,int m) {
  double d = 0.0;
  for (int j=0; j<m; j++)
    d += square(w[j]-v[j]);
  return sqrt(d);
}
#endif
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static int dofmap_mult(Mat A,Vec x,Vec y) { // y =A*x
  void *ctx;
  int ierr = MatShellGetContext(A,&ctx); 
  CHKERRQ(ierr); 
  Dofmap * dofmap = (Dofmap *)ctx;
  dofmap->mult(A,x,y);
  return 0;
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int Dofmap::mult(Mat A,Vec x,Vec y) { // y =A*x
  int ierr;
  int nrow = nnod*ndof;
  int ncol = neqtot;
  double *xp,*yp;
  ierr = VecGetArray(x,&xp);
  assert(!ierr);
  ierr = VecGetArray(y,&yp);
  assert(!ierr);

  w.resize(nrow);
  
  for (int j=0; j<nrow; j++) w[j]=0.0;
  qxpy(&w[0],xp,1.0);
  for (int j=0; j<neqtot; j++) yp[j]=0.0;
  qtxpy(yp,&w[0],1.0);

  ierr = VecRestoreArray(x,&xp);
  assert(!ierr);
  ierr = VecRestoreArray(y,&yp);
  assert(!ierr);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Dofmap::solve(double *xp,double *yp) {
  int nrow = nnod*ndof;
  int ncol = neqtot;
  int m, ierr;

  // We have to solve the system Q*x = y. 
  // In this version we solve it with the CG
  // on the normal equations. Normally the matrix Q
  // is a permutation and in other cases it is very well
  // conditioned so that this is a good choice. So we
  // solve Q'*Q*x = Q'*y, -> H*x = z

  // Store x, z in PETSc vectors
  ierr = VecCreateSeq(PETSC_COMM_SELF,
		      neqtot,&x); 
  assert(!ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF,
		      neqtot,&z); 
  assert(!ierr);

  // `w' is an auxiliary vector, contains `Q*x'. 
  w.resize(nrow);

  // Compute `z = Q*y'
  double *zp;
  ierr = VecGetArray(z,&zp); 
  assert(!ierr);
  qtxpy(zp,yp,1.0);
  ierr = VecRestoreArray(z,&zp);
  assert(!ierr);

  // Define auxiliary matrix shell
  Mat A;
  ierr = MatCreateShell(PETSC_COMM_SELF,ncol,
			ncol,ncol,ncol,this,&A);
  assert(!ierr);
  MatShellSetOperation(A,MATOP_MULT,
		       (void (*)(void))(&dofmap_mult));
  assert(!ierr);

  // Define auxiliary KSP
  KSP ksp;
  PC pc;
  ierr = KSPCreate(PETSC_COMM_SELF,&ksp); assert(!ierr);
  ierr = KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN); assert(!ierr);
  ierr = KSPGetPC(ksp,&pc); assert(!ierr);
  ierr = KSPSetType(ksp,KSPCG); assert(!ierr);
  ierr = PCSetType(pc,PCNONE); assert(!ierr);
  // ierr = KSPSetMonitor(ksp,KSPDefaultMonitor,NULL,NULL);

  // Solve problem
  int its;
  ierr = KSPSolve(ksp,z,x);
  ierr = KSPGetIterationNumber(ksp,&its); assert(!ierr);

  assert(its<=100);
  printf("Dofmap::solve: solved projection in %d iters\n",its);

#if 0
  double a=1.0;
  ierr = VecSet(&a,x);  assert(!ierr);
  ierr = MatMult(A,x,z);  assert(!ierr);
  ierr = VecView(z,PETSC_VIEWER_STDOUT_SELF);
  assert(!ierr);
#endif

  // Copy auxiliary PETSc vector `x' in
  // output argument `x'
  double *xpp;
  ierr = VecGetArray(x,&xpp); 
  assert(!ierr);
  for (int j=0; j<neqtot; j++)
    xp[j] = xpp[j];
  ierr = VecRestoreArray(x,&xpp);
  assert(!ierr);

  // Destroy auxiliary quantities
  ierr =MatDestroy(&A); assert(!ierr);
  ierr = VecDestroy(&x); assert(!ierr);
  ierr = VecDestroy(&z); assert(!ierr);
  ierr = KSPDestroy(&ksp); assert(!ierr);
  w.clear();
}
