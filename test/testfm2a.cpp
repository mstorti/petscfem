/*__INSERT_LICENSE__*/
//$Id: testfm2a.cpp,v 1.3 2001/05/30 03:58:53 mstorti Exp $
 
#include <stdio.h>
#include "../src/fastmat2.h"

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

int main() {

  FastMat2 A(2,3,3);
  init123(A);
  A.print("A: ");
  A.ir(1,1);
  Indx indx;
  A.get_dims(indx);
  indx.print("free indices");
  A.print("A first row: ");

  A.rs().print("A: " );
  A.ir(2,2).print("A column 2: ");
  A.rs().is(2,3).print("Column 3");
  A.is(2).print("all the matrix: ");
  A.rs();

  A.ir(2,3).print("Column 3 (but as row now)");
  A.ir(2).print("all the matrix: ");
  A.rs();

  A.resize(2,5,5).set(0.);
  FastMat2 B;
  // set some elements and print
  A.setel(1.,1,1);
  A.setel(2.,1,2);
  A.setel(3.,2,1);
  A.setel(4.,2,2);
  A.print("A: ");

  // B is the trace of A
  B.sum(A,-1,-1);
  double a = double(B);

  // Adds to a element
  A.addel(1.,1,1);
  A.addel(2.,1,2);
  A.addel(3.,2,1);
  A.addel(4.,2,2);
  A.print("Adds to elements and gets 2*A: ");

  // Divides by 2.
  A.scale(0.5).print("Scales by 0.5 and gets the same back:");

  // select a subblock
  A.r(2).r(1).c(2).c(1).print("Sub-block (2-1)(2-1)");

  // transpose
  A.t().print("A: ");

  // copy to B
  B.resize(2,2,2).set(A).print("B: ");

  // make row, print, and back to 2x2 matrix
  B.reshape(1,4).print("B (reshaped): ");
  B.reshape(2,2,2);

  // transpose back
  A.t().print("A': ");

  // set to B transpose
  B.t().set(A).print("B': ");
  
  // add A to B
  B.add(A).print("B+A (2*A): ");

  // rest A to B, so that it remains the same
  B.rest(A).print("B+A-A (B): ");

  // add 10*A to B
  B.axpy(A,10.).print("B+10.*A: ");

  // set B to a and summ 2*A
  B.set(4.).axpy(A,2.).print("4+2.*A");

  // reset A, set column 1 to 2
  A.rs().set(0.).c(1).set(2.).print("A set to 2. (column1) ");

  FastMat2 AA(2,2,3),BB(2,3,2),C,CC;
  AA.set(1.).print("A: ");
  BB.set(2.).print("B: ");
  C.prod(AA,BB,1,-1,-1,2).print("C = A*B : ");
  AA.ir(1,1);
  BB.ir(1,1);
  CC.prod(AA,BB,1,2).print("Kron prod of AA(1,:),BB(1,:)");

  FastMat2 D(4,3,2,2,3),E,F;
  D.set(4.);
  D.print("D: ");
  E.ctr(D,-1,1,2,-1).print("D_{ijki}: ");
  F.prod(D,D,-1,-2,-3,1,-1,-2,-3,2).print("F_{lm} = D_{ijkl} D_{ijkm}");
  E.print("E = D_{ijjk} = ");

  FastMat2 G(3,2,2,2),H;
  G.set(-3.).print("G: ");
  H.sum(G,-1,1,2).print("H = Sum over first index of G: ");
  H.sum_square(G,-1,1,2).print("H = Sum_square over first index of G: ");
  H.sum_abs(G,-1,1,2).print("H = Sum_abs over columns of G: ");

  G.setel(-5.,1,1,1).setel(5.,2,2,2).print("G: ");
  H.min(G,-1,1,2).print("H = Min over first index of G: ");
  H.max(G,-1,1,2).print("H = Max over first index of G: ");
  H.min_abs(G,-1,1,2).print("H = Min abs over first index of G: ");
  H.max_abs(G,-1,1,2).print("H = Max abs over first index of G: ");

  FastMat2 K(3,2,3,3),M,N;
  init123(K);
  K.print("K: ");
  M.diag(K,1,-2,-2).print("M = diag(K) : ");

}
