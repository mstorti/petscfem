//__INSERT_LICENSE__
//$Id: testfm2.cpp,v 1.4 2001/04/01 01:35:07 mstorti Exp $

#include <stdio.h>
#include <time.h>

#include <sles.h>

#include "../src/fastmat2.h"
#include "../src/fem.h"
#include "../src/utils.h"
#include "../src/util2.h"
#include <newmatio.h>

void init123(FastMat2 & A) {
  Indx dims,indx;
  A.get_dims(dims);
  int ndims = dims.size();
  indx = Indx(ndims,1);
  while (1) {
    double val=0.;
    for (int jj=0; jj<ndims; jj++) 
      val = val*10+indx[jj];
    A.setel(indx,val);
    if (!inc(indx,dims)) break;
  }
}

double dexp (double x) {
  double y = exp(x);
  return y;
}

struct mydata_t {
  double lambda; 
  int n;
};
  

double myfun(double x,void *user_args) {

  mydata_t *mydata_p = (mydata_t *) user_args;
  return pow((x / mydata_p->lambda),double(mydata_p->n));

}

int main() {

  Chrono chrono;
  FastMatCacheList cache_list;
  FastMatCachePosition cp1,cp2;
  //int n=3,Nin=1000,Nout=10;
  int n=3,Nin=10,Nout=10;
  FastMat2 AA(2,n,n),A(2,3,3),B,C,D,E(2,3,3),G,H,K,L,M,P,Q,R,S,T,U,V,
    W,X,XX(2,n,n),Y,Z,Z1,Z2(1,3),Z3,Z4(2,3,3),Z121,Z20(2,3,5),Z21(3,3,5,5);
  FastMat2 Z5(2,20,20),Z6;
  FMatrix Z7(3,3),Z8,Z9,Z10(2,2),Z11,Z12,Z15(4,4),Z16,Z17,Z18,
    Z116,Z117,Z118,Z19(3,3);
  Matrix NA(3,3),NB;
  NA << 1. << 3. << 5. << 7. << 9. << 11. << 13. << 15. << 17;
  A.set(NA);
  init123(AA);
  init123(Z5);
  init123(XX);
  init123(Z15);
  Z15.setel(20.,1,1).setel(20.,2,2);
  double e[] = {1.,2.,4.,5.,3.,2.,3.,2.,1.};
  double pp[9],dz8,z12,z121,z13,z14,z15,z115;
  chrono.start();
  srand(time(0));
  double z10[] = {1.,2.,3.,4.};
  Z10.set(z10);
  init123(Z20);
  init123(Z21);
  
  mydata_t mydata;
  mydata.lambda=10.;
  mydata.n=2;
  
  for (int j=0; j<Nout; j++) {
    FastMat2::activate_cache(&cache_list);
    for (int k=0; k<Nin; k++) {
      FastMat2::reset_cache();

      A.set(AA);
      B.set(A);

      // Setting individual elements
      A.setel(1.,1,1);
      A.setel(2.,1,2);
      A.setel(3.,2,1);
      A.setel(4.,2,2);

      // element-by-element multiplication (Matlab's .*)
      B.mult(A);

#define CP(cp)	printf("\n==========\n ------> Enter code portion [ " \
		       #cp " ]\n==========\n\n")

      FastMat2::branch();
      if (drand()<.8) {

	//	CP(0);
	FastMat2::choose(0);

	FastMat2::branch();
	if (drand()<.5) {

	  //	  CP(0.0);
	  FastMat2::choose(0);
	  // Copying matrices
	  C.set(A);

	  // Adding to individual addition
	  C.addel(.5,1,1);
	  C.addel(.5,1,2);
	  C.addel(.5,2,1);
	  C.addel(.5,2,2);

	  // Max (as generic sum) of all elements
	  z14 = A.min_all();

	} else {

	  //	  CP(0.1);
	  FastMat2::choose(1);
	  D.set(A);

	  // Multiplying to individual addition
	  D.multel(.5,1,1);
	  D.multel(.5,1,2);
	  D.multel(.5,2,1);
	  D.multel(.5,2,2);

	  // Max (as generic sum) of all elements
	  z13 = A.max_all();

	}
	FastMat2::leave();

	// Setting from an array of doubles
	E.set(e);

	// Addition on an index (Matlab's sum())
	G.sum(A,1,-1);

	// Sum of all elements
	z12 = A.sum_square_all();
	z121 = Z121.sum_square(A,-1,-1).get();


      } else {
	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

	//	CP(1);
	FastMat2::choose(1);
	// Sum of squares and abs() over a column
	H.sum_square(A,1,-1);
	K.set(A).setel(-1000.,2,2);
	L.sum_abs(K,1,-1);


	// Several maximum 
	M.max(K,-1,1);
	P.min(K,-1,1);
	Q.max_abs(K,-1,1);
	R.min_abs(K,-1,1);
      
	// Set all elements
	S.set(A).set(3.14);

	// Scale elements
	U.set(A).scale(2.);

	// Add to elements
	V.set(A).add(5.);

	// Product of matrices.
	W.prod(A,A,1,-1,-1,2);

	// Reshaping 
	XX.reshape(1,n*n);
	X.set(XX);
	XX.reshape(2,n,n);

	// Importing from and exporting to Newmat
	Y.set(NA);
	Y.export_vals(NB);
      
      // Taking the diagonal
	Z1.diag(A,-1,-1);

      }
      FastMat2::leave();
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

      // getting individual elements.
      double a,b,c;
      a=A.get(1,1);
      b=A.get(2,1);
      c=A.get(1,2);
      Z2.setel(a,1).setel(b,2).setel(c,3);

      // Contraction (trace)
      Z3.ctr(A,-1,-1);

      // Testing a loop (a cache is created for each execution of the loop)
      for (int j=0; j<3; j++) {
	Z4.ir(2,3-j).set(A.ir(2,j+1));
      }
      A.rs();
      Z4.rs();

      // Testing indices with more than 10 elements
      Z5.rs().is(1,1,5).is(2,2,7);
      Z5.rs().is(1,1,5).is(1,2,20).is(2,3,4);
      // Z6.sum(Z5,-1,-1);

      Z7.set(e);
      Z8.set(Z7);

      // Test det and inverse
      dz8 = Z8.det();
      Z9.inv(Z8);

      // Kronecker 
      Z11.kron(A,Z10);

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // The following loop is executed a varying number of times
      // but it is executed at least one time. 
      FastMat2::get_cache_position(cp1);
      int max_loop=5;
      int n=irand(1,max_loop);
      for (int ll=0; ll<n; ll++) {
	FastMat2::jump_to(cp1);
	// The following operations are done indirectly via Newmat Matrices
	// Determinant for dim>3
	z15 = Z15.det();

	// Inverse for dim>3
	Z16.inv(Z15);

	// Apply a scalar function element by element
	Z17.set(Z16).fun(&dexp);
      }
      FastMat2::resync_was_cached();

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // Apply a scalar function element by element with auxiliary arguments
      Z18.set(Z16).fun(&myfun,(void *)&mydata);

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // Now it is executed may be zero times. This poses a problem if
      // the first time it is executed zero times.
      FastMat2::branch();
      FastMat2::choose(0);

      FastMat2::get_cache_position(cp2);
      n=irand(0,max_loop-1);
      if (k==0) n=0; // This is the critical case. n=0 the first execution
		     // of the loop. 
      for (int ll=0; ll<n; ll++) {
	FastMat2::jump_to(cp2);
	// The following operations are done indirectly via Newmat Matrices
	// Determinant for dim>3
	z115 = Z15.det();

	// Inverse for dim>3
	Z116.inv(Z15);

	// Apply a scalar function element by element
	Z117.set(Z116).fun(&dexp);
      }
      FastMat2::resync_was_cached();
      FastMat2::leave();

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // Apply a scalar function element by element with auxiliary arguments
      Z118.set(Z16).fun(&myfun,(void *)&mydata);

      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // Set to identity
      Z19.set(5.34).eye(3.342);

      Z21.d(2,3).set(Z20).rs();
    }
    FastMat2::void_cache();
  }
  FastMat2::deactivate_cache();
  double cpu=chrono.elapsed();
  printf("copy (with set): n=%d, Nin,Nout=%d,%d, elapsed: %f, "
	 "speed: %f Mflops \n",
	 n,Nin,Nout,cpu,FastMat2::operation_count()*double(Nin*Nout)/1e6/cpu);
  
#define SH(n) n.print(#n ": ")

  SH(A);
  SH(B);
  SH(C);
  SH(D);
  SH(E);
  SH(G);
  SH(K);
  SH(L);
  SH(M);
  SH(P);
  SH(Q);
  SH(R);
  SH(S);
  SH(U);
  SH(V);
  SH(W);
  SH(X);
  SH(Y);
  SHV(NA);
  SHV(NB);
  SH(Z);
  SH(Z1);
  SH(Z2);
  SH(Z3);
  SHV(double(Z3));
  SH(Z5);
  SH(Z6);
  SH(Z8);
  SHV(dz8);
  SH(Z9);
  SH(Z11);
  SH(Z12);
  SHV(z12);
  SHV(z121);
  SHV(z13);
  SHV(z14);
  SHV(z15);
  SH(Z16);
  SH(Z17);
  SH(Z18);
  SHV(z115);
  SH(Z116);
  SH(Z117);
  SH(Z118);
  SH(Z19);
  Z21.d(3,2);SH(Z21);
  Z21.rs().d(2,3);SH(Z21);

#undef SH

}
