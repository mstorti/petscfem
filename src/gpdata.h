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
 
#ifndef GPDATA_H
#define GPDATA_H
#include <stdio.h>
#include <newmatio.h>

#include "fastmat.h"
#include "fastmat2.h"

// Options for GPdata
#define GP_NEWMAT 0
#define GP_FASTMAT 1
#define GP_FASTMAT2 2

/** Contains information on Gauss integration points (GP). 
    @author M. Storti
*/ 
class GPdata {
public:
  /** gradient of shape functions with respect to master element
      coordiantes at each GP
  */
  Matrix *dshapexi,
    /** gradient of shape functions with respect to global
	coordiantes at each GP
    */
    *dshapex;
  /// For the FastMat version
  FastMat **FM_dshapexi, **FM_shape;
  /// For the FastMat2 version
  FastMat2 **FM2_dshapexi, **FM2_shape;
  /// Shape functions
  RowVector *shape;
  /// Integrations weights
  double *wpg;
  /** Constructor. 
      @author M. Storti
      @param geom (input) string defining the geometry, currently may be
      cartesian<n>d, with n=1,2,3, triangle or tetra. 
      @param ndim (input) number of dimensions. (This may be
      different from the problem dimension. Typically it may be one
      dimension lower for surface elements.)
      @param nel (input) number of nodes per element
      @param npg (input) number of gauss points
      @param mat_version (input) 
  */ 
   GPdata(const char *geom,int ndim,int nel,int npg,int mat_version=GP_NEWMAT);
  /// Destructor
  ~GPdata (void);
  /// number of Gauss points. 
  int npg;
  /// flags whether Newmat or FastMat
  int mat_version;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Provides the shape functions and gradients for the bilinear
    quadrangle. 
    @author M. Storti
    @param shape (output) the shape function of the different nodes at
    the Gauss Points. 
    @param dshapexi (output) the gradients of the shape functions with
    respect to master element coordinates. 
    @param xipg (input) xi coordinate of the Gauss point. 
    @param etapg (input) eta coordinate of the Gauss point. 
*/ 
int cartesian_2d_shape(RowVector &shape,Matrix &dshapexi,
		       double xipg,double etapg);

#endif
