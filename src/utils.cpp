//__INSERT_LICENSE__
//$Id: utils.cpp,v 1.9 2001/08/02 01:54:01 mstorti Exp $
 
#include <stdio.h>

#include <set>
#include <cassert>

#include <newmatio.h>

#include "fem.h"
#include "utils.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "mydet"
double mydet(Matrix A) {
  int m = A.Nrows();
  if (m==1) {
    return A(1,1);
   } else if(m==2) {
     return A(1,1)*A(2,2)-A(1,2)*A(2,1);
   } else if (m==3) {
     return A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
       + A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))
       + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1));
   } else {
     printf("Not implemented order (%d)  for determinant at file: %s, line: %d\n",
	    m,__FILE__,__LINE__);
   }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "mydetsur"
/** Computes the area of the hipersurface defined by the columns of
    Jaco. Jaco is the jacobian from ndimel intrinsic coordinates in
    the master element to the ndim cartesian coordinates, so that 
    the j-th row of Jaco represents \dep{!x}{\xi_j} = e_j. We should
    compute the generalized vectorial product of the rows of jaco and
    compute the length of this resulting vector. 
 */
double mydetsur(Matrix & Jaco, ColumnVector &S) {
  int ierr;
  int m=Jaco.Nrows(), n=Jaco.Ncols();
  assert(m==n-1);
  // if (m!=n-1) PFEMERRQ("m!=n-1. (m=%d, n=%d)\n");
  Matrix JJ(n,n);
  //  RowVector S(n);
  //  S.ReSize(n);
  JJ.SubMatrix(1,m,1,n) = Jaco;
  for (int k=1; k<=n; k++) {
    JJ.SubMatrix(n,n,1,n) = 0;
    JJ(n,k) = 1;
    S(k) = mydet(JJ);
  }
  // fixme:= I put this in order to fix a bug but it should be coded
  // more cleanly if it will be used in the future. 
  S=-S; 
  double SS = sqrt(S.SumSquare());
  return SS;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "double mydetsur(FastMat2 &Jaco, FastMat2 &S)"
// fixme:= I put this in order to fix a bug but it should be coded
// more cleanly if it will be used in the future. 
#if 1 
double mydetsur(FastMat2 &Jaco, FastMat2 &S) {
  int ierr;
  int n=Jaco.dim(2);
  assert(Jaco.dim(1)==n-1);
  double sx,sy,sz;

  if (n==2) {
    Jaco.ir(1,1);
    sx=Jaco.get(2);
    sy=-Jaco.get(1);
    S.setel(sx,1).setel(sy,2);
    Jaco.rs();
    return sqrt(sx*sx+sy*sy);
  } else if (n==3) {
    sx = -Jaco.get(1,2)*Jaco.get(2,3)+Jaco.get(1,3)*Jaco.get(2,2);
    sy = -Jaco.get(1,3)*Jaco.get(2,1)+Jaco.get(1,1)*Jaco.get(2,3);
    sz = -Jaco.get(1,1)*Jaco.get(2,2)+Jaco.get(1,2)*Jaco.get(2,1);

    S.setel(sx,1).setel(sy,2).setel(sz,3);
    return sqrt(sx*sx+sy*sy+sz*sz);
  } else {
    printf("Dimension not allowed (only 2 and 3)");
    assert(0);
  }
}
#else  // old version
double mydetsur(FastMat2 &Jaco, FastMat2 &S) {
  int ierr;
  int n=Jaco.dim(2);
  assert(Jaco.dim(1)==n-1);
  double sx,sy,sz;

  if (n==2) {
    Jaco.ir(1,1);
    sx=-Jaco.get(2);
    sy=Jaco.get(1);
    S.setel(sx,1).setel(sy,2);
    Jaco.rs();
    return sqrt(sx*sx+sy*sy);
  } else if (n==3) {
    sx = Jaco.get(1,2)*Jaco.get(2,3)-Jaco.get(1,3)*Jaco.get(2,2);
    sy = Jaco.get(1,3)*Jaco.get(2,1)-Jaco.get(1,1)*Jaco.get(2,3);
    sz = Jaco.get(1,1)*Jaco.get(2,2)-Jaco.get(1,2)*Jaco.get(2,1);

    S.setel(sx,1).setel(sy,2).setel(sz,3);
    return sqrt(sx*sx+sy*sy+sz*sz);
  } else {
    printf("Dimension not allowed (only 2 and 3)");
    assert(0);
  }
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "kron"
Matrix kron(const Matrix & A,const Matrix & B) {
  int na,ma,nb,mb,nc,mc;
  ma=A.Nrows();
  na=A.Ncols();
  mb=B.Nrows();
  nb=B.Ncols();
  mc=ma*mb;
  nc=na*nb;

  Matrix C(mc,nc);

  for (int j=0; j<ma; j++) {
    for (int k=0; k<na; k++) {
      int j1,j2,k1,k2;
      j1=j*mb+1;
      j2=j1+mb-1;
      k1=k*nb+1;
      k2=k1+nb-1;
      C.SubMatrix(j1,j2,k1,k2) = A(j+1,k+1)*B;
    }
  }
  return C;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0 // Already defined as inline in utils.h
double drand() {
  return ((double)(rand()))/((double)(RAND_MAX));
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int mini(int n,...) {
  va_list list;
  va_start(list,n);
  int min,item;
  for (int kk=0; kk<n; kk++) {
    item = va_arg(list,int);
    min = ( kk==0 ? item : ( min < item ? min : item));
  }
  va_end(list);
  return min;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int maxi(int n,...) {
  va_list list;
  va_start(list,n);
  int max,item;
  for (int kk=0; kk<n; kk++) {
    item = va_arg(list,int);
    max = ( kk==0 ? item : ( max > item ? max : item));
  }
  va_end(list);
  return max;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double maxd(int n,...) {
  va_list list;
  va_start(list,n);
  double max,item;
  for (int kk=0; kk<n; kk++) {
    item = va_arg(list,double);
    max = ( kk==0 ? item : ( max > item ? max : item));
  }
  va_end(list);
  return max;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "reshape(Matrix &A,int m,int n)"
int reshape(Matrix &A,int m,int n) {
  int mm=A.Nrows();
  int nn=A.Ncols();
  PFEMERRCQ(mm*nn != m*n,"Total number of elements don't coincide");
  Matrix B(m,n);
  B << A.Store();
  A = B;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int wait_from_console(char *s=NULL) {
  static int deac=0;
  if (deac) return 0;
  int myrank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);
  int ierr;
  char ans;
  ierr = MPI_Barrier(PETSC_COMM_WORLD);
  if (myrank==0) {
    if (s!=NULL) printf("%s --- ",s);
    printf("Continue? (n/RET=y) > ");
    fflush(stdout);
    scanf("%c",&ans);
  }
  ierr = MPI_Bcast (&ans, 1, MPI_CHAR, 0,PETSC_COMM_WORLD);
  CHKERRQ(ierr); 
  if (ans=='n') {
    PetscFinalize();
    exit(0);
  } else if (ans=='d') {
    deac = 1;
  } 
  ierr = MPI_Barrier(PETSC_COMM_WORLD);
  CHKERRQ(ierr);  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Estos templates no se si van en los header o en los .cpp (???)
#if 0
template<class T> 
class T random_pop(set<T> &Tset) {
  T retval;
  int n=Tset.size();
  int j = irand(1,n);
  set<T>::iterator it=Tset.begin();
  for (int k=1; k<j; k++) it++;
  retval = *it;
  Tset.erase(it);
  return retval;
}
#endif
