// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: gpdata.h,v 1.11 2004/01/26 20:22:34 mstorti Exp $
 
#ifndef GPDATA_H
#define GPDATA_H
#include <stdio.h>
#include <newmatio.h>

#include <src/fastmat.h>
#include <src/fastmat2.h>
#ifdef USE_DX
#include <src/util3.h>
#endif

// Options for GPdata
#define GP_NEWMAT 0
#define GP_FASTMAT 1
#define GP_FASTMAT2 2

/** Contains information on Gauss integration points (GP). 
    @author M. Storti
*/ 
class GPdata {
public:
  /** Gradient of shape functions with respect to master element
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
  /// The volume of the reference (master) element
  double master_volume;
  /// Sum of Integrations weights
  double wpg_sum;
#ifdef USE_DX
  /** This represents how the element is mapped onto Data
      Explorer connections. */
  DXSplit splitting;
#endif

  friend class edge;

private:
  int nedges_m;
  vector<int> edges;
public:
  class edge {
  private:
    friend class GPdata;
    int indx;
    const GPdata *gp;
    edge(int j,const GPdata *gp_a);
  public:
    edge();
    edge operator++(int);
    friend bool operator!=(edge p, edge q);
    int first();
    int second();
  };
  edge edges_begin();
  edge edges_end();
  int nedges();

public:
  static int get_default_geom(int ndim, int nel, string& geometry);
  static int get_default_npg(const string& geometry, int &npg);

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
