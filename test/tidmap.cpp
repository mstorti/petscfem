/*
  This file belongs to the PETSc - FEM package, a library and
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
 
#include <newmatio.h>
#include "../src/utils.h"
#include "../src/idmap.h"

// Tests the `idmap' class

void scale(row_t & row, const double c) {
  row_t::iterator it;
  for (it=row.begin(); it!=row.end(); it++) {
    if (c==0.) {
      row.erase(it);
    } else {
      it->second *= c;
    }
  }
}

void check (int n,idmap &id, Matrix &ID, double &maxerr) {
  maxerr=0.;
  for (int j=1; j<=n; j++) {
    for (int k=1; k<=n; k++) {
      double err;
      id.get_val(j,k,err);
      err -= ID(j,k);
      err=fabs(err);
      if (err>maxerr) maxerr=err;
    }
  }
}

int check_idmap(int n, int N, int verbose_print) {

  printf("testing for size of idmap n=%d, number of operations %d\n",n,N);
  double tol=1e-10;
  idmap id(n,IDENTITY_MAP);
  if (verbose_print) id.print("before permutations\n");
  row_t r1,r2;

  Matrix ID(n,n);
  RowVector R1,R2;
  ColumnVector C1,C2;
  for (int j=1; j<=n; j++)
    ID(j,j)=1.;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Permutation of rows check
  for (int j=0; j<N; j++) {
    int i1=irand(1,n);
    int i2=irand(1,n);
    
    printf("exchanging rows %d, %d\n",i1,i2);
    
    // Intercambia filas en idmap
    id.get_row(i1,r1);
    id.get_row(i2,r2);
    id.del_row(i1);
    id.del_row(i2);
    id.row_set(i1,r2);
    id.row_set(i2,r1);

    // intercambia en ID
    R1=ID.Row(i1);
    R2=ID.Row(i2);
    ID.Row(i1)=R2;
    ID.Row(i2)=R1;

  }
  if (verbose_print) {
    id.print("after row permutations.\n");
    cout << ID;
  }

  // checks internal consistency
  id.check();

  // Checks that they are equal
  double maxerr=0.;
  for (int j=1; j<=n; j++) {
    for (int k=1; k<=n; k++) {
      double err;
      id.get_val(j,k,err);
      err -= ID(j,k);
      err=fabs(err);
      if (err>maxerr) maxerr=err;
    }
  }
  printf("Row permutations: maximum error: %f \nOK? : %s\n",maxerr,
	 (maxerr<tol ? "YES" : "NOT"));

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Permutation of cols check
  for (int j=0; j<N; j++) {
    int i1=irand(1,n);
    int i2=irand(1,n);
    
    printf("exchanging cols %d, %d\n",i1,i2);
    
    // Intercambia filas en idmap
    id.get_col(i1,r1);
    id.get_col(i2,r2);
    id.del_col(i1);
    id.del_col(i2);
    id.column_set(i1,r2);
    id.column_set(i2,r1);

    // intercambia en ID
    C1=ID.Column(i1);
    C2=ID.Column(i2);
    ID.Column(i1)=C2;
    ID.Column(i2)=C1;

  }
  if (verbose_print) {
    id.print("After column permutations.\n");
    cout << ID;
  }

  // checks internal consistency
  id.check();

  // Checks that they are equal
  maxerr=0.;
  for (int j=1; j<=n; j++) {
    for (int k=1; k<=n; k++) {
      double err;
      id.get_val(j,k,err);
      err -= ID(j,k);
      err=fabs(err);
      if (err>maxerr) maxerr=err;
    }
  }
  printf("Column permutations: maximum error: %f, \nOK? : %s\n",maxerr,
	 (maxerr<tol ? "YES" : "NOT"));

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Linear combination of rows
  row_t r;
  for (int j=0; j<N; j++) {
    int i1,i2;
    while (1) {
      i1=irand(1,n);
      i2=irand(1,n);
      if (i1!=i2) break;
    }
    double c1=2*drand()-1.;
    double c2=2*drand()-1.;
    double c = sqrt(c1*c1+c2*c2);
    c1 /= c;
    c2 /= c;
    
    printf("rotating rows %d, %d, %f, %f\n",i1,i2,c1,c2);
    
    // Intercambia filas en idmap
    id.get_row(i1,r1);
    id.get_row(i2,r2);
    id.del_row(i1);
    id.del_row(i2);

    VOID_IT(r);
    r = r1;
    scale(r,c1);
    axpy(r,c2,r2);
    id.row_set(i1,r);

    VOID_IT(r);
    r = r2;
    scale(r,c1);
    axpy(r,-c2,r1);
    id.row_set(i2,r);

    // intercambia en ID
    R1=ID.Row(i1);
    R2=ID.Row(i2);
    ID.Row(i1)= c1*R1 + c2*R2;
    ID.Row(i2)=-c2*R1 + c1*R2;

  }
  if (verbose_print) {
    id.print("after row rotations.\n");
    cout << ID;
  }

  // checks internal consistency
  id.check();

  // Checks that they are equal
  maxerr=0.;
  for (int j=1; j<=n; j++) {
    for (int k=1; k<=n; k++) {
      double err;
      id.get_val(j,k,err);
      err -= ID(j,k);
      err=fabs(err);
      if (err>maxerr) maxerr=err;
    }
  }
  printf("Row rotations: maximum error: %f \nOK? : %s\n",maxerr,
	 (maxerr<tol ? "YES" : "NOT"));

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Linear combination of columns
  for (int j=0; j<N; j++) {
    int i1,i2;
    while (1) {
      i1=irand(1,n);
      i2=irand(1,n);
      if (i1!=i2) break;
    }
    double c1=2*drand()-1.;
    double c2=2*drand()-1.;
    double c = sqrt(c1*c1+c2*c2);
    c1 /= c;
    c2 /= c;
    
    printf("rotating columnss %d, %d, %f, %f\n",i1,i2,c1,c2);
    
    // Intercambia filas en idmap
    id.get_col(i1,r1);
    id.get_col(i2,r2);
    id.del_col(i1);
    id.del_col(i2);

    VOID_IT(r);
    r = r1;
    scale(r,c1);
    axpy(r,c2,r2);
    id.column_set(i1,r);

    VOID_IT(r);
    r = r2;
    scale(r,c1);
    axpy(r,-c2,r1);
    id.column_set(i2,r);

    // intercambia en ID
    C1=ID.Column(i1);
    C2=ID.Column(i2);
    ID.Column(i1)= c1*C1 + c2*C2;
    ID.Column(i2)=-c2*C1 + c1*C2;

  }
  if (verbose_print) {
    id.print("after column rotations.\n");
    cout << ID;
  }

  // checks internal consistency
  id.check();

  // Checks that they are equal
  check (n,id,ID,maxerr);
  printf("Column rotations: maximum error: %f \nOK? : %s\n",maxerr,
	 (maxerr<tol ? "YES" : "NOT"));

  // Testing set_elem
  for (int k=0; k<N; k++) {
    int i=irand(1,n);
    int j=irand(1,n);
    double val=2*drand()-1.;
    printf("setting (%d,%d) to %f\n",i,j,val);
    id.set_elem(i,j,val);
    ID(i,j)=val;
  }

  check (n,id,ID,maxerr);
  printf("After set_elem(): maximum error: %f \nOK? : %s\n",maxerr,
	 (maxerr<tol ? "YES" : "NOT"));
  
  // Testing get_val()
  maxerr=0;
  for (int k=0; k<N; k++) {
    int i=irand(1,n);
    int j=irand(1,n);
    double val,err;
    id.get_val(i,j,val);
    printf("testing (%d,%d): from idmap=%f, from ID=%f\n",i,j,val,ID(i,j));
    err = fabs(val-ID(i,j));
    if (err>maxerr) maxerr=err;
  }
  printf("After get_val(): maximum error: %f \nOK? : %s\n",maxerr,
	 (maxerr<tol ? "YES" : "NOT"));

  // Tests solve()
  for (int j=1; j<=n; j++) {
    id.del_row(j);
  }
  ID=0;
  
  for (int j=1; j<=n; j++) {
    id.set_elem(j,j,1.);
    ID(j,j)=1.;
  }

  double *x = new double[n];
  double *y = new double[n];
  ColumnVector Y(n),X(n);
  for (int k=1; k<n/3; k++) {
    int i=irand(1,n);
    int j=irand(1,n);
    double c=2.*drand()-1.;
    id.set_elem(i,j,c);
    ID(i,j)=c;
  }

  for (int k=1; k<=n; k++) {
    y[k-1] = 2.*drand()-1.;
    Y(k) = y[k-1];
  }

  for (int j=0; j<N; j++) {
    int i1=irand(1,n);
    int i2=irand(1,n);
    
    printf("exchanging rows %d, %d\n",i1,i2);
    
    // Intercambia filas en idmap
    id.get_row(i1,r1);
    id.get_row(i2,r2);
    id.del_row(i1);
    id.del_row(i2);
    id.row_set(i1,r2);
    id.row_set(i2,r1);

    // intercambia en ID
    R1=ID.Row(i1);
    R2=ID.Row(i2);
    ID.Row(i1)=R2;
    ID.Row(i2)=R1;

  }

  if(verbose_print) {
    id.print("matrix for solve()\n");
    cout << ID;
  }

  check (n,id,ID,maxerr);
  printf("After permuting and setting random: maximum error: %f \nOK? : %s\n",maxerr,
	 (maxerr<tol ? "YES" : "NOT"));

  id.solve(x,y);
  X = ID.i()*Y;

  maxerr=0;
  for (int j=1; j<=n; j++) {
    if (verbose_print) printf("%d -> %f %f \n",j,X(j),x[j-1]);
    double err= fabs(X(j)-x[j-1]);
    if (err>maxerr) maxerr=err;
  }
  printf("After matrix solving: maximum error: %f \nOK? : %s\n",maxerr,
	 (maxerr<tol ? "YES" : "NOT"));

  delete[] x;
  delete[] y;

}

int main() {

  check_idmap(5,100,1);
  check_idmap(100,50,0);

}
