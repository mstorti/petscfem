/* $Id: walldata.cpp,v 1.1.1.1 2000/12/28 12:54:43 mstorti Exp $ */

/*
  This file belongs to he PETSc - FEM package a library and
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
 
#include "../../src/fem.h"
//  #include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/utils.h"
#include "../../src/util2.h"
#include "nsi_tet.h"

extern int MY_RANK,SIZE;


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "WallData::WallData()"
WallData::WallData(vector<double> *data_pts_,vector<ElemToPtr>
		   *elemset_pointer_,int ndim_) {
  ndim=ndim_;
  npoints = data_pts_->size()/ndim;
  data_pts = annAllocPts(npoints,ndim);

  for (int k=0; k<npoints; k++) {
    for (int j=0; j<ndim; j++) {
      data_pts[k][j] = (*data_pts_)[k*ndim+j];
    }
  }
  data_pts_->resize(0);
  kd_tree = new ANNkd_tree(data_pts,npoints,ndim);

  // copies temporary vector<ElemToPtr> to ElemToPtr[]
  // This is tricky?
  nelemset = elemset_pointer_->size();
  elemset_pointer = elemset_pointer_->begin();
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void WallData::nearest()"
void WallData::nearest(const ANNpoint &point, Elemset *& elemset, int &elem, ANNidx &nn_idx,
			ANNpoint &nn,ANNdist &dist) {
  static int KNBR=1;
  kd_tree->annkSearch(point,KNBR,&nn_idx,&dist,0);
  for (int k=0; k<ndim; k++) 
    nn[k] = data_pts[nn_idx][k];

  int kk;
  for (kk=0; kk<nelemset; kk++) {
    if (nn_idx< elemset_pointer[kk].first) break;
  }
  elemset = elemset_pointer[kk].second;
  int prev=0;
  if (kk>0) prev = elemset_pointer[kk-1].first;
  elem = nn_idx-prev;
}
    


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void WallData::nearest()"
void WallData::nearest(double *point,int &nn) {

  static int KNBR=1;
  double dist;
  kd_tree->annkSearch(point,KNBR,&nn,&dist,0);
}
    

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void WallData::nearest_elem_info"
void WallData::nearest_elem_info(const int nn, Elemset *& elemset,
				 int &elem, const double *& coords) {
  int kk;
  for (kk=0; kk<nelemset; kk++) {
    if (nn< elemset_pointer[kk].first) break;
  }
  elemset = elemset_pointer[kk].second;
  int prev=0;
  if (kk>0) prev = elemset_pointer[kk-1].first;
  elem = nn-prev;
  coords = data_pts[kk];
}
    
