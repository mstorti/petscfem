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
 
#include <iostream>
#include <stdio.h>
#include "../src/fastlib.h"

int main() {
  // dimension a vector of 10 components
  int n=10;
  FastVector<int> a(n,0);
  for (int j=0; j<n; j++) {
    a[j] = j;
  }

  // resizes to 40
  int nn=40;
  a.resize(nn);
  for (int j=0; j<nn; j++) {
    a[j] = j;
  }
  a.print("resized");

  // resizes back to 20
  nn=20;
  a.resize(nn);
  for (int j=0; j<nn; j++) {
    a[j] = j;
  }
  a.print("resized");

  FastVector<double> c(n,3.3);
  c.print("with doubles");

}
