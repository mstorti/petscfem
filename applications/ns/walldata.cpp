//__INSERT_LICENSE__
//$Id: walldata.cpp,v 1.13 2006/02/18 22:40:47 mstorti Exp $
 
#include <src/fem.h>
//  #include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/sttfilter.h>
#include "nsi_tet.h"

extern int MY_RANK,SIZE;


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "WallData::WallData()"
WallData::WallData() 
  : kd_tree(NULL), data_pts(NULL) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void WallData
::init(vector<double> *data_pts_a,
       vector<ElemToPtr> *elemset_pointer_a,int ndim_) {
  clear();
#ifdef USE_ANN
  ndim=ndim_;

  /// The position of the points
  npoints = data_pts_a->size()/ndim;
  data_pts = annAllocPts(npoints,ndim);

  for (int k=0; k<npoints; k++) {
    for (int j=0; j<ndim; j++) {
      data_pts[k][j] = (*data_pts_a)[k*ndim+j];
    }
  }
  // This should go to wall::after_assemble()
  data_pts_a->resize(0);
  kd_tree = new ANNkd_tree(data_pts,npoints,ndim);

  // copies temporary vector<ElemToPtr> to ElemToPtr[]
  // This is tricky?
  nelemset = elemset_pointer_a->size();
  elemset_pointer = &*elemset_pointer_a->begin();
#else
  PETSCFEM_ERROR0("Not compiled with ANN library!!\n");
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void WallData::clear() {
  if (kd_tree) delete kd_tree;
  kd_tree = NULL;
  // fixme:= hmmmmm.... who deletes data_pts???
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "WallData::WallData()"
WallData::~WallData() { clear(); }


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#ifdef USE_ANN
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
#endif


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void WallData::nearest()"
void WallData::nearest(double *point,int &nn) {

#ifdef USE_ANN
  static int KNBR=1;
  double dist;
  kd_tree->annkSearch(point,KNBR,&nn,&dist,0);
#else
  PETSCFEM_ERROR0("Not compiled with ANN library!!\n");
#endif
}
    

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void WallData::nearest_elem_info"
void WallData::nearest_elem_info(const int nn, Elemset *& elemset,
				 int &elem, const double *& coords) {
#ifdef USE_ANN
  int kk;
  for (kk=0; kk<nelemset; kk++) {
    if (nn< elemset_pointer[kk].first) break;
  }
  elemset = elemset_pointer[kk].second;
  int prev=0;
  if (kk>0) prev = elemset_pointer[kk-1].first;
  elem = nn-prev;
  // coords = data_pts[kk]; // wrong!!
  coords = data_pts[nn];
#else
    PETSCFEM_ERROR0("Not compiled with ANN library!!\n");
#endif
}
