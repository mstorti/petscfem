// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: nsi_tet.h,v 1.12 2001/06/29 20:19:22 mstorti Exp $
#ifndef NSI_TET_H  
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
//  class nsi_tet : public ns_volume_element { 
//  public: 
//    ASSEMBLE_FUNCTION;
//  };

//  //-------<*>-------<*>-------<*>-------<*>-------<*>------- 
//  class nsi_tet_les : public ns_volume_element { 
//  public: 
//    ASSEMBLE_FUNCTION;
//  };

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_fm2 : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_les_ther : public ns_volume_element { 
public: 
  ASSEMBLE_FUNCTION;
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
class nsi_tet_keps : public ns_volume_element { 
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
class bcconv_nsther_fm2 : public Elemset { 
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Global parameters used by NS
struct GlobParam {
  /// Trapezoidal rule integration parameter
  double alpha;
  /// time step
  double Dt;
  /// do not include temporal term (steady solutions)
  int steady;
};

void wall_fun(double yp,double &f,double &fprime);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class WallFun {
public:
  virtual void init() {};
  virtual void w(double yp,double &f,double &fp)=0;
  virtual ~WallFun()=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class WallFun1 : public WallFun {
  Elemset *elemset;
public:
  void w(double yp,double &f,double &fprime);
  // assuming a wall law of the form f = 2.5*log(yplus) + 5.5
  // then if it has to be compatible with f = 1/Chi * log(E*yplus)
  // we find: Chi = 0.4; E = 9.025;
  // WallFun1(double Chi_=0.4, double E_=9.025) : Chi(Chi_), E_star(E_) {};
  // WallFun1() : Chi(Chi_), E_star(E_) {};
  WallFun1(Elemset *e) : elemset(e) {};
  ~WallFun1() {};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class NonLinearRes : public Elemset {
 public:
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
  virtual int nres()=0;
  virtual void init()=0;
  virtual void res(int k,FastMat2 &U,FastMat2 & r,
		   FastMat2 & lambda,FastMat2 & jac)=0;
  virtual ~NonLinearRes()=0;
};

#define MAXPROP_WLR 10
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class wall_law_res : public NonLinearRes {
  int nk,ne,ndim;
  double y_wall_plus,fwall,fprime,C_mu,
    viscosity,von_Karman_cnst,coef_k,coef_e,
    turbulence_coef;
  WallFun *wf;
  int u_wall_indx,nprops;
  FastMat2 u_wall,du_wall;
  int elprpsindx[MAXPROP_WLR]; 
  double propel[MAXPROP_WLR];
public:
  int nres() {return 2;};
  void init();
  void res(int k,FastMat2 &U,FastMat2 & r,FastMat2 & lambda,
	   FastMat2 & jac);
  wall_law_res() {wf = new WallFun1(this);};
  ~wall_law_res() {delete wf;};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class wallke : public Elemset { 
public: 
  ASK_FUNCTION;
  ASSEMBLE_FUNCTION;
};

#endif
