#ifndef NSI_TET_H  // -*- mode: C++ -*-
#define NSI_TET_H

#include <ANN/ANN.h>			// ANN declarations
#include <vector>			// ANN declarations

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// This is the typical element for volume computations. The `ask'
// function is the same.
class ns_volume_element : public Elemset { 
public: 
  ASK_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_fm2 : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class bcconv_ns_fm2 : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class wall : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

typedef pair<int,Elemset *> ElemToPtr;

typedef vector<ElemToPtr> ElemToPtrV;

class WallData {
private:
  /// The octree
  ANNkd_tree *kd_tree;                 // search structure
  /// The position of the points
  ANNpointArray	data_pts;		// data points
  /// Number of points
  int npoints;
  /// number of spacial dimensions
  int ndim;
  /// pointers to the corresponding elemsets
  ElemToPtr *elemset_pointer;
  /// The length of the elemset_pointer array
  int nelemset;
public:
  /// constructor from a vector of coordinates, pointers and dimensions
  WallData(vector<double> *data_pts_,vector<ElemToPtr>
	   *elemset_pointer,int ndim_);
  /// find the nearest neighbor
  void nearest(const ANNpoint &point, Elemset *& elemset, int &elem, ANNidx &nn_idx,
	       ANNpoint &nn,ANNdist &dist);
  /// find the nearest neighbor (short version) only returns the index.
  void nearest(double *point,int &nn);
  /** Give wall element info. 
      Given the nearest neighbor index give the corresponding
      wall element info: elemset, element number, coords of
      center of wall element. 
      @author M. Storti
      @param nn (input) the index of the wall element in
      the WallData structure as returned by nearest()
      @param elemset (output) the `wall' elemset to which the element pertains
      @param elem (output) the wall element index (base 0)
      @param coords (output) the coords of the center of the wall element. 
  */
  void WallData::nearest_elem_info(const int nn, Elemset *& elemset, int &elem,
				   const double *& coords);
};

#endif

