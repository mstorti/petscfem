// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: surf2vol2.h,v 1.1 2003/02/25 03:14:04 mstorti Exp $
#ifndef PETSCFEM_SURF2VOL2_H
#define PETSCFEM_SURF2VOL2_H

#include <src/surf2vol.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Instantiation of the #Surf2Vol# class for
    hexas with quad surface elements. */
class Quad2Hexa : public Surf2Vol {
private:
  /** Stores the three possible combinations of faces
      (not taking into account reflections about
      mid-planes and rotations abut axes. */
  static const int faces[][8];
  /// Auxiliary array for element connectivity
  int vol[8];
  /// Auxiliary array for element connectivity
  int vol_r[8];
  /// Auxiliary array for surface connectivity
  int this_face[4];
public:
  /// Constructor. Sets #use_exterior_normal_m#
  Quad2Hexa(const char *geom,int ndim,int nel,
	   int npg,int mat_version=GP_NEWMAT,
	    int use_exterior_normal_m=0) 
    : Surf2Vol(geom,ndim,nel,npg,mat_version,
	       use_exterior_normal_m) { 
    assert(!use_exterior_normal_m);
  }
  /** @name Callback routines for the quad/hexa combination. */
  //@{
  /** Callback routine that defines the possible
      orientations of a face. 
      @param j (input) number of face
      @param fc (output) connectivity of the #j#-th surface face. 
      @param vol (input) rotation map corresponding to the
      #j#-th volume */
  void face(int j,const int *&fc,const int *&vol);
  /** Returns number of elements, surface nodes and volume nodes. */
  int nfaces(int &nel_surf,int &nel_vol) { nel_surf=4; nel_vol=8; return 24; }
  // void surface_nodes(int &nel_surf,int &nel_vol);
  //@}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Instantiation of the #Surf2Vol# class for
    prisms with tri surface elements. */
class Tri2Prism : public Surf2Vol {
private:
  int vol_c[6],vol[6];
  void rotate(int *fc,int nrot);
  void reflect(int *fc);
  void invert(int *fc);
public:
  /// Constructor. Sets #use_exterior_normal_m#
  Tri2Prism(const char *geom,int ndim,int nel,
	   int npg,int mat_version=GP_NEWMAT,
	    int use_exterior_normal_m=0) 
    : Surf2Vol(geom,ndim,nel,npg,mat_version,
	       use_exterior_normal_m) { }
  /** @name Callback routines for the quad/hexa combination. */
  //@{
  /** Callback routine that defines the possible
      orientations of a face. 
      @param j (input) number of face
      @param fc (output) connectivity of the #j#-th surface face. 
      @param vol (input) rotation map corresponding to the
      #j#-th volume */
  void face(int j,const int *&fc,const int *&vol);
  /** Returns number of elements, surface nodes and volume nodes. */
  int nfaces(int &nel_surf,int &nel_vol) { nel_surf=3; nel_vol=6; return 6; }
  // void surface_nodes(int &nel_surf,int &nel_vol);
  //@}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Instantiation of the #Surf2Vol# class for
    quads with line surface elements. */
class Line2Quad : public Surf2Vol {
public:
  /// Constructor. Sets #use_exterior_normal_m#
  Line2Quad(const char *geom,int ndim,int nel,
	   int npg,int mat_version=GP_NEWMAT,
	    int use_exterior_normal_m=0) 
    : Surf2Vol(geom,ndim,nel,npg,mat_version,
	       use_exterior_normal_m) { }
  /** @name Callback routines for the line/quad combination. */
  //@{
  /** Callback routine that defines the possible
      orientations of a face. 
      @param j (input) number of face
      @param fc (output) connectivity of the #j#-th surface face. 
      @param vol (input) rotation map corresponding to the
      #j#-th volume */
  void face(int j,const int *&fc,const int *&vol);
  /** Returns number of elements, surface nodes and volume nodes. */
  int nfaces(int &nel_surf,int &nel_vol) { nel_surf=2; nel_vol=4; return 4; }
  // void surface_nodes(int &nel_surf,int &nel_vol);
  //@}
};

#endif
