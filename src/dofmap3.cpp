//__INSERT_LICENSE__
//$Id: dofmap3.cpp,v 1.7 2005/02/20 15:37:31 mstorti Exp $

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
  int ncol = neq;
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
  int ncol = neq;
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Dofmap::solve(double *x,double *y) {
  int nrow = nnod*ndof;
  int ncol = neq;
  int m;
  vector<double> zv(nrow,0.0), 
    vv(nrow,0.0), wv(ncol,0.0),
    dv(ncol,0.0);
  double *z = &zv[0], *v = &vv[0],
    *w = &wv[0], *d = &dv[0];
  int niter=100;
  double rtol=1e-10,atol=1e-10,r0;
  double omega=1;
  // Make product z=Qy
  qtxpy(z,y,1.0);
  // printv(z,nrow,"z");
  // assert(nrow==ncol);
  // printf("|z-y| %f\n",diff(z,y,nrow));
  // Compute diagonal
  for (int j=0; j<nrow; j++) {
    int kdof = j % ndof + 1;
    int node = (j / ndof) + 1;
    const int *dofs;
    const double *coef;
    get_row(node,kdof,m,&dofs,&coef);
    for (int k=0; k<m; k++) {
      int dof = dofs[k]-1;
      d[dof] += square(coef[k]);
    }
#if 0
    if (m>1 || (dofs[0]-1)!=j)
      for (int k=0; k<m; k++)
	printf(" %d %d  %f\n",j,dofs[k]-1,coef[k]);
#endif
  }
  // printv(d,nrow,"d");
  for (int iter=0; iter<niter; iter++) { 
    // v = 0.0;
    for (int j=0; j<nrow; j++) v[j]=0.0;
    // v = Q*x
    qxpy(v,x,1.0);
    // printf("|v-x| %f\n",diff(v,x,nrow));
    // printv(v,ncol,"v");
    // w = Q'*Q*x
    for (int j=0; j<ncol; j++) w[j]=0.0;
    qtxpy(w,v,1.0);
    // printf("|w-v| %f\n",diff(w,v,nrow));
    // printv(w,nrow,"w");
    // printf("|z-w| %f\n",diff(z,w,nrow));
    double res=0.0;
    for (int j=0; j<ncol; j++) {
      x[j] += omega*(z[j]-w[j])/d[j];
      res += square(z[j]-w[j]);
    }
    // printv(x,ncol,"x");
    res = sqrt(res);
    if (iter==0) r0 = res;
    if (res<atol || res<rtol*r0) break;
    printf("iter %d, res %f\n",iter,res);
  }
}
