// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: surf2vol.h,v 1.3 2003/02/25 20:34:22 mstorti Exp $
#ifndef PETSCFEM_SURF2VOL_H
#define PETSCFEM_SURF2VOL_H

#include <src/gpdata.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This is a Gauss Point data class but for elements that match a
    surface element with a volume element.  The typical case is for
    integrating a surface term that depends on normal gradients of the
    state variables. For instance, when computing viscous forces. The
    user enters only the connectivities of the surface elements and
    the `initialize()' function finds the volume element that shares
    a face with this element. */
class Surf2Vol : public GPdata {
private:
  /// Flags whether to use the exterior normal or the interior normal
  int use_exterior_normal_m;
public:
  /** Constructor as for the `GPdata' class plus the flag
      `use_exterior_normal'. 
      @param geom (input) the geometry of the integration element
      @param ndim (input) the dimension of the volume element. 
      @param nel (input) number of elements in the volume element. 
      @param mat_version (input) the matrix version (Newmat or FastMat2)
      @param use_exterior_normal (input) flags whether we the faces are
      oriented so that the normal to it points towards the interior or
      the exterior of the exterior of the volume element. */
  Surf2Vol(const char *geom,int ndim,int nel,int npg,
	   int mat_version=GP_NEWMAT,int use_exterior_normal=0);

  /** @name Call back functions. */ 
  //@{
  /** For each interpolation family we have to provide the number of
      possible faces (including rotations orientations). For instance
      for an hexahedral element we have 6x4=24 possible faces.  
      @param nel_surf (input) number of nodes in the surface element
      @param nel_vol (input) number of nodes in the volume element
      @return  the number of possible faces. */
  virtual int nfaces(int &nel_surf,int &nel_vol)=0;
  /** Returns the $j$-th possible face. 
      @param j (input) the of face number (0 based)
      @param fc (output) an array of #nel_surf# nodes (0 based)
      that represent the connectivity of the $j$-th face
      @param vol (input) an array of #nel_vol# integers that
      represent the mapping of the connectivity of the volume
      element in order to rotate it to a standard position. */ 
  virtual void face(int j,const int *&fc,const int *&vol)=0;
  //@}

  /** Rotates the connectivity in #vol_map# according to
      the surface rotation #surf_map#.
      @param surf_map (input) the mapping of the surface
      nodes to the local nodes in the voume element (both 0 based). 
      @param vol_map (input/output) remaps the connectivity of
      volume element (#nel_vol# nodes) so that it remains in
      standard position. */
  int map_mask(const int *surf_map,int *vol_map);
  /** Callback ruotines may ask with this routines 
      @return #true# if callback routines must return face
      connectivities so that their normal is oriented to
      the exterior normal and #false# otherwise. */
  int use_exterior_normal() { return (use_exterior_normal_m); }

  static void factory(TextHashTable *thash, string &volume_elemset,
		      Surf2Vol *sv_gp_data, Elemset *& vol_elem,
		      int &identify_volume_elements,int &layers,
		      int &use_exterior_normal,int &ndimel);

};

void identify_volume_elements_fun(int nnod, int nel_surf, int layers,
				  int nelem, int *icone, int nel, int nel_vol,
				  Elemset *vol_elem,Surf2Vol *sv_gp_data);

#endif
