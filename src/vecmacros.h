// -*- mode: c++ -*-
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

#ifndef VECMACROS_H
#define VECMACROS_H

/** @name Vec macros to allow Fortran-like access to array elements. */
//@{

/// adressing a 1-dimensional array as 2 dimensional. 
#define VEC_ADDR_2(j,k,dk) (dk)*(j)+(k)
/// adressing a 1-dimensional array as 3 dimensional. 
#define VEC_ADDR_3(j,k,dk,l,dl) (VEC_ADDR_2(j,k,dk))*(dl)+(l)
/// adressing a 1-dimensional array as 4 dimensional. 
#define VEC_ADDR_4(j,k,dk,l,dl,p,dp) (VEC_ADDR_3(j,k,dk,l,dl))*(dp)+(p)
/// adressing a 1-dimensional array as 5 dimensional. 
#define VEC_ADDR_5(j,k,dk,l,dl,p,dp,q,dq) \
                       (VEC_ADDR_4(j,k,dk,l,dl,p,dp))*(dq)+(q)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Accessing a 1-dimensional array as 2 dimensional. 
    Typical use: \#define MATRIX(j,k) VEC2(j,k,dk)
    @author M. Storti
    @param name (input) name of the array
    @param j (input) row index 
    @param k (input) column index (running faster)
    @param dk (input) column dimension
*/ 
#define VEC2(name,j,k,dk) ((name)[VEC_ADDR_2(j,k,dk)])

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Accessing a 1-dimensional array as 3 dimensional. 
    Last index runs faster.
    Typical use: \#define MATRIX(j,k,l) VEC2(j,k,dk,l,dl)
    @author M. Storti
    @param name (input) name of the array
    @param j (input) first index 
    @param k (input) 2nd index
    @param dk (input) 2nd index dimension
    @param l (input) 3rd index
    @param dl (input) 3rd index dimension
*/ 
#define VEC3(name,j,k,dk,l,dl) ((name)[VEC_ADDR_3(j,k,dk,l,dl)])

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Accessing a 1-dimensional array as 4 dimensional. 
    Last index runs faster.
    Typical use: \#define MATRIX(j,k,l) VEC2(j,k,dk,l,dl)
    @author M. Storti
    @param name (input) name of the array
    @param j (input) first index 
    @param k (input) 2nd index
    @param dk (input) 2nd index dimension
    @param l (input) 3rd index
    @param dl (input) 3rd index dimension
    @param p (input) 4th index
    @param dp (input) 4th index dimension
*/ 
#define VEC4(name,j,k,dk,l,dl,p,dp) ((name)[VEC_ADDR_4(j,k,dk,l,dl,p,dp)])

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Accessing a 1-dimensional array as 4 dimensional. 
    Last index runs faster.
    Typical use: \#define MATRIX(j,k,l,p,q) VEC2(j,k,dk,l,dl,p,dp,q,dq)
    @author M. Storti
    @param name (input) name of the array
    @param j (input) 1st index 
    @param k (input) 2nd index
    @param dk (input) 2nd index dimension
    @param l (input) 3rd index
    @param dl (input) 3rd index dimension
    @param p (input) 4th index
    @param dp (input) 4th index dimension
    @param q (input) 5th index
    @param dq (input) 5th index dimension
*/ 
#define VEC5(name,j,k,dk,l,dl,p,dp,q,dq) \
           ((name)[VEC_ADDR_5(j,k,dk,l,dl,p,dp,q,dq)])
//@}
#endif
