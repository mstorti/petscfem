//__INSERT_LICENSE__
//$Id: fastmat.cpp,v 1.5 2001/05/30 18:21:53 mstorti Exp $

#include <math.h>
#include <stdio.h>
#include <cassert>

#include <newmat.h>
#include "fastmat.h"

int FastMat::count=0;

//#define  DEBUG_FASTMAT_ALLOC

#ifdef FASTMAT_DEBUG_DESTR
int FastMat::lastid=0;
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

FastMat::FastMat(int m_=0,int n_=0,double *s=NULL) :
  m(m_), n(n_), store(s) {

  m=0;
  n=0;
  defined=0;
  FM_COUNT(++count);
#ifdef FASTMAT_DEBUG_DESTR
  ident=0;
#endif

  if (m_==0 || n_==0) return;
  set_size(m_,n_);
  if (s!=NULL) set(s);
#if 0
  defined=1;
#ifdef FASTMAT_DEBUG_DESTR
  ident=++lastid;
#endif
  store = new double[m*n];
#ifdef DEBUG_FASTMAT_ALLOC
  printf("FastMat alloc %d %d\n",m,n);
#endif
  if (s==NULL) return;
  double *to=store;
  double *from=s;
  for (int k=0;k<m*n;k++) {
    *(to++) = *(from++);
  }
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::set_size(int m_,int n_) {
  assert(m==0 && n==0 && store==NULL && !defined);
  m=m_;
  n=n_;
#ifdef DEBUG_FASTMAT_ALLOC
  printf("FastMat alloc %d %d\n",m,n);
#endif
  store = new double[m*n];
  defined=1;
#ifdef FASTMAT_DEBUG_DESTR
  ident=++lastid;
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::print(char *s=NULL) {
  assert(defined);
  printf("%s\n",(s==NULL ? "values: " : s));
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      printf("%g ",*(store+i*n+j));
    }
    printf("\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::set(int i,int j,double val=0) {
#ifdef FASTMAT_BASE_1
  *(store+(i-1)*n+j-1) = val;
#else
  *(store+i*n+j) = val;
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:     
/// adds val to element i,j
void FastMat::add(const int i,const int j,const double val) {
  *(LOCATION(i,j)) += val;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::set(double *val,int nn=0) {
  if (nn==0) nn=m*n;
  double *from=val,*to=store;
  for (int i=0; i<nn; i++) {
    *(to++) = *(from++);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::set(const FastMat & A,int nn=0) {
  if (!defined) set_size(A.m,A.n);
  set(A.store,nn);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::copy(double *val,int nn=0) const {
  if (nn==0) nn=m*n;
  double *from=store,*to=val;
  for (int i=0; i<nn; i++) {
    *(to++) = *(from++);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::set(double val) {
  int nn=m*n;
  double *to=store;
  for (int i=0; i<nn; i++) {
    *(to++) = val;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::get(const int i,const int j,double * val) const {
#ifdef FASTMAT_BASE_1 
  *val = *(store+(i-1)*n+j-1);
#else
  *val = *(store+i*n+j);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::get(const int i,const int j,double & val) const {
#ifdef FASTMAT_BASE_1
  val = *(store+(i-1)*n+j-1);
#else
  val = *(store+i*n+j);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat::~FastMat() {
  assert(defined || store==NULL);
#ifdef DEBUG_FASTMAT_ALLOC
  if (defined) printf("FastMat delete %d %d\n",m,n);
#endif
  if (defined) delete[] store;
  FM_COUNT(--count);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat::eye(int m_=0) {

  // it is defined or m_ non null
  assert(m>0 || m_>0);

  // check if not defined
  if (!defined) set_size(m_,m_);

  // check dimensions
  assert(m==m_ && n==m_);
  
  set(0.);
  double *to = store;
  int inc=m+1;
  for (int i=0; i<m; i++) {
    *to = 1.;
    to += inc;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FMp(FastMat & A,const FastMat & B,const FastMat & C) {

  // Set size (If not set already)
  if (A.m==0) A.set_size(B.m,C.n);
  
  // Check matrix sizes
  assert(B.n == C.m && A.m == B.m && A.n == C.n);

  for (int i=0; i<A.m; i++) {
    for (int j=0; j<A.n; j++) {
      double val=0;
      double *b,*c;
      b = B.store+i*B.n;
      c = C.store+j;
      for (int k=0; k<B.n; k++) {
	val += (*b) * (*c);
	b++;
	c += C.n;
      }
      *(A.store+i*A.n+j) = val;
    }
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Matrix Addition
int FMa(FastMat & A,const FastMat & B,const FastMat & C) {

  // Set size (If not set already)
  if (A.m==0) A.set_size(B.m,B.n);
  
  assert(B.n == C.n && A.n == C.n);
  assert(B.m == C.m && A.m == C.m);

  double *a=A.store;
  double *b=B.store;
  double *c=C.store;

  for (int i=0; i<(A.m*A.n); i++) {
    *(a++) = *(b++) + *(c++);
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FMaxpy(FastMat & Y,double a,FastMat & X) {
  assert(X.m == Y.m && X.n == Y.n );

  double *x=X.store;
  double *y=Y.store;

  for (int i=0; i<(X.m*X.n); i++) {
    *(y++) += a * *(x++);
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FMaxpy_t(FastMat & Y,double a,const FastMat & X) {
  assert(X.m == Y.n && X.n == Y.m );

  double *from,*to;
  for (int i=0; i<Y.m; i++) {
    from = X.location0(0,i);
    to = Y.location0(i,0);
    for (int j=0; j<Y.n; j++) {
      *to += a*(*from);
      to++;
      from += X.n;
    }
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FMinv(FastMat & invA,FastMat & A) {
#define AA(i,j) VEC2(A.store,i,j,m)
#define IAA(i,j) VEC2(invA.store,i,j,m)

  // check matrix is square
  assert(A.m==A.n);

  // Set size (If not set already)
  if (invA.m==0) invA.set_size(A.m,A.m);

  // check dimensions
  assert(invA.m == A.m && invA.n == A.m);
  
  int m = A.m;
  double det;
  int ierr = FMdet(det,A); FMCHK;
  if (det==0) {
    printf("Determinant is null.\n");
    return 1;
  }
  if (m==1) {
    invA.store[0] = 1./A.store[0];
    return 0;
  } else if (m==2) {
    IAA(0,0) = AA(1,1)/det;
    IAA(0,1) = -AA(0,1)/det;
    IAA(1,0) = -AA(1,0)/det;
    IAA(1,1) = AA(0,0)/det;
    return 0;
  } else if (m==3) {
    IAA(0,0) = (AA(1,1)*AA(2,2)-AA(1,2)*AA(2,1))/det;
    IAA(0,1) = (AA(2,1)*AA(0,2)-AA(2,2)*AA(0,1))/det;
    IAA(0,2) = (AA(0,1)*AA(1,2)-AA(0,2)*AA(1,1))/det;
	   
    IAA(1,0) = (AA(1,2)*AA(2,0)-AA(1,0)*AA(2,2))/det;
    IAA(1,1) = (AA(2,2)*AA(0,0)-AA(2,0)*AA(0,2))/det;
    IAA(1,2) = (AA(0,2)*AA(1,0)-AA(0,0)*AA(1,2))/det;
	   
    IAA(2,0) = (AA(1,0)*AA(2,1)-AA(1,1)*AA(2,0))/det;
    IAA(2,1) = (AA(2,0)*AA(0,1)-AA(2,1)*AA(0,0))/det;
    IAA(2,2) = (AA(0,0)*AA(1,1)-AA(0,1)*AA(1,0))/det;
    return 0;
  } else {
    printf("Not allowed dimension.\n");
    return 1;
  }
}
    
#undef AA
#undef IAA
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FMtr(FastMat & At,const FastMat & A) {
  double *a,*at;
  if (At.n != A.m || At.m != A.n ) return 1;
  for (int i=0; i<A.m; i++) {
    a=A.store+i*A.n;
    at=At.store+i;
    for (int j=0; j<A.n; j++) {
      //      printf("At(%d,%d) -> %f\n",j,i,*at);
      *at = *(a++);
      at += A.m;
    }
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FMdet(double & det,const FastMat & A) {
  if (A.m != A.n) return 1;
  int m =A.n;
  double *store = A.store;
#define AA(i,j) VEC2(store,i,j,m)
   if (m==1) {
     det = AA(1,1);
     return 0;
   } else if(m==2) {
     det = AA(0,0)*AA(1,1)-AA(0,1)*AA(1,0);
     return 0;
   } else if (m==3) {
     det = AA(0,0)*(AA(1,1)*AA(2,2)-AA(1,2)*AA(2,1))
       + AA(0,1)*(AA(1,2)*AA(2,0)-AA(1,0)*AA(2,2))
       + AA(0,2)*(AA(1,0)*AA(2,1)-AA(1,1)*AA(2,0));
     return 0;
   } else {
     printf("Not implemented order (%d)  for determinant.");
     return 1;
   }
#undef AA
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
FastMat::FastMat(const Matrix & B) {

#ifdef FASTMAT_DEBUG_DESTR
  ident=++lastid;
#endif
  m=B.Nrows();
  n=B.Ncols();
  defined=1;
  set(B.Store());
  print();
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int kron(FastMat & C, FastMat const & A, FastMat const & B) {
  
  // set if defined
  if (C.m==0) C.set_size(A.m*B.m,A.n*B.n);

  // check dimensions
  assert(C.m == A.m*B.m);
  assert(C.n == A.n*B.n);
  
  for (int i=0; i<A.m; i++) {
    for (int j=0; j<A.n; j++) {
      double a;
      A.get(i,j,a);
      for (int k=0; k<B.m; k++) {
	double *from,*to;
	from = B.store+k*B.n;
	to = C.store+(((i*B.m+k)*A.n+j)*B.n);
	for (int l=1; l<=B.n; l++) *(to++) = *(from++)*a;
      }
    }
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Get some columns
void FastMat::get(int i1,int i2,int j1,int j2,const FastMat & B) {

  // dimension if not already set
  if (!defined) {
    assert(j2>=j1);
    assert(i2>=i1);
    set_size(i2-i1+1,j2-j1+1);
  }
  
  // Check dimensions
  assert(m==i2-i1+1);
  assert(n==j2-j1+1);
  
  double *to,*from;
  for (int i=0; i<m; i++) {
    to = location0(i,0);
    from = B.location0(i1+i-FASTMAT_BASE,j1-FASTMAT_BASE);
    for (int j=0; j<n; j++) *(to++) = *(from++);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Get some columns
void FastMat::columns(const FastMat & B,int j1,int j2) {

  // dimension if not already set
  if (!defined) {
    assert(j2>=j1);
    set_size(B.m,j2-j1+1);
  }
  
  // Check dimensions
  assert(m==B.m);
  assert(n==j2-j1+1);
  
  double *to,*from;
  for (int i=0; i<m; i++) {
    to = location0(i,0);
    from = B.location0(i,j1-FASTMAT_BASE);
    for (int j=0; j<n; j++) *(to++) = *(from++);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Get one column
void FastMat::column(const FastMat & B,int j) {

  // dimension if not already set
  if (!defined) {
    assert(1<=j && j<=B.n);
    set_size(B.m,1);
  }
  
  // Check dimensions
  assert(m==B.m);
  assert(n==1);
  
  double *to,*from;
  to = location0(0,0);
  from = B.location0(0,j-FASTMAT_BASE);
  for (int i=0; i<m; i++) {
    *(to++) = *from;
    from += B.n;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Get some rows
void FastMat::rows(const FastMat & B,int i1,int i2) {

  // dimension if not already set
  if (!defined) {
    assert(i2>=i1);
    set_size(i2-i1+1,B.n);
  }
  
  // Check dimensions
  assert(m==i2-i1+1);
  assert(n==B.n);
  
  double *to,*from;
  for (int i=0; i<m; i++) {
    to = location0(i,0);
    from = B.location0(i1+i-FASTMAT_BASE,0);
    for (int j=0; j<n; j++) *(to++) = *(from++);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Get one row
void FastMat::row(const FastMat & B,int i) {

  // dimension if not already set
  if (!defined) {
    assert(1<=i && i<=B.m);
    set_size(1,B.n);
  }
  
  // Check dimensions
  assert(m==1);
  assert(n==B.n);
  
  double *to,*from;
  to = location0(0,0);
  from = B.location0(i-FASTMAT_BASE,0);
  for (int j=0; j<n; j++) *(to++) = *(from++);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// convert from Newmat matrices 
int NM2FM(FastMat & A,const Matrix & B) {
  if (!A.is_defined()) A.set_size(B.Nrows(),B.Ncols());
  A.set(B.Store());
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// reshapes the matrix (keeping size)
void FastMat::reshape(int m_, int n_) {
  if (m==0) set_size(m_,n_);
  assert(m*n==m_*n_);
  m=m_;
  n=n_;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// sum of squares
int FastMat::sum_square(double & val) const {
  double c,sum=0;
  double * from = store;
  int ii=m*n;
  for (int i=1; i<=ii; i++) {
    c = *(from++);
    sum += c*c;
  }
  val = sum;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// sum of squares
double FastMat::sum() const {
  double c,val=0;
  double * from = store;
  int ii=m*n;
  for (int i=1; i<=ii; i++) {
    c = *(from++);
    val += c;
  }
  return val;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FastMat::transpose(const FastMat & B) {

  if (!defined) set_size(B.n,B.m);
  assert(m==B.n && n==B.m);

  double *from,*to;
  for (int i=0; i<m; i++) {
    from=B.location0(0,i);
    to=location0(i,0);
    int inc=B.n;
    for (int j=0; j<n; j++) {
      *to = *from;
      to++;
      from += inc;
    }
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FastMat::scale(const double c) {

  double * from = store;
  int ii=m*n;
  for (int i=1; i<=ii; i++) *(from++) *= c;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int FastMat::trace_of_product(const FastMat & B,double & trace) const {
  double s=0;
  assert(n==B.m && m==B.n);
  double *fromA=store,*fromB=B.store;
  int nn=B.n;
  for (int i=0; i<m; i++) {
    fromA = location0(i,0);
    fromB = B.location0(0,i);
    for (int j=0; j<n; j++) {
      s += (*fromA) * (*fromB);
      fromA++;
      fromB += nn;
    }
  }
  trace = s;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "double FastMat::norm1() const" 
double FastMat::norm1() const {
  double max=0,sum,*from;
  int inc=n;
  for (int j=0; j<n; j++) {
    sum=0;
    from = location0(0,j);
    for (int i=0; i<m; i++) {
      sum += fabs(*from);
      from += inc;
    }
    if (sum>max) max=sum;
  }
  return max;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// sets submatrix  i1 to i2, j1 to j1 to val
void FastMat::set(int i1,int i2, int j1, int j2,double *val) {

  assert(defined && i2>=i1 && j2>=j1 && i2<m+FASTMAT_BASE  &&
	 j2<n+FASTMAT_BASE);
  
  double *from=val, *to;
  for (int i=i1; i<=i2; i++) {
    to = LOCATION(i,j1);
    for (int j=j1; j<=j2; j++) {
      *(to++) = *(from++);
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// sets submatrix  i1 to i2, j1 to j1 to matrix B
void FastMat::set(int i1,int i2, int j1, int j2,const FastMat & B) {
  assert(defined && i2>=i1 && j2>=j1 && i2<m+FASTMAT_BASE  &&
	 j2<n+FASTMAT_BASE);
  assert(i2-i1+1 == B.m && j2-j1+1 == B.n);
  set(i1,i2,j1,j2,B.store);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// sets submatrix  i1 to i2, j1 to j1 to val
void FastMat::add(int i1,int i2, int j1, int j2,double *val) {

  assert(defined && i2>=i1 && j2>=j1 && i2<m+FASTMAT_BASE  &&
	 j2<n+FASTMAT_BASE);
  
  double *from=val, *to;
  for (int i=i1; i<=i2; i++) {
    to = LOCATION(i,j1);
    for (int j=j1; j<=j2; j++) {
      *(to++) += *(from++);
    }
  }
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// sets submatrix  i1 to i2, j1 to j1 to val
void FastMat::add(int i1,int i2, int j1, int j2,const FastMat & B) {

  assert(defined && i2>=i1 && j2>=j1 && i2<m+FASTMAT_BASE  &&
	 j2<n+FASTMAT_BASE);
  assert(i2-i1+1 == B.m && j2-j1+1 == B.n);
  add(i1,i2,j1,j2,B.store);
#if 0  
  double *from=B.store, *to;
  for (int i=i1; i<=i2; i++) {
    to = LOCATION(i1,j1);
    for (int j=j1; j<=j2; j++) {
      *(to++) += *(from++);
    }
  }
#endif
}
