// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: embgath.h,v 1.13 2003/01/06 03:07:43 mstorti Exp $
#ifndef EMBGATH_H
#define EMBGATH_H

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
      the exterior of the exterior of the volume element. 
  */
  Surf2Vol(const char *geom,int ndim,int nel,int npg,
	   int mat_version=GP_NEWMAT,int use_exterior_normal=0);

  /** @name Call back functions. 
  */ 
  //@{
  /** For each interpolation family we have to provide the number of
      possible faces (including rotations orientations). For instance
      for an hexahedral element we have 6x4=24 possible faces.  
      @param nel_surf (input) number of nodes in the surface element
      @param nel_vol (input) number of nodes in the volume element
      @return  the number of possible faces.
  */
  virtual int nfaces(int &nel_surf,int &nel_vol)=0;
  /** Returns the $j$-th possible face. 
      @param j (input) the of face number (0 based)
      @param fc (output) an array of #nel_surf# nodes (0 based)
      that represent the connectivity of the $j$-th face
      @param vol (input) an array of #nel_vol# integers that
      represent the mapping of the connectivity of the volume
      element in order to rotate it to a standard position. 
  */ 
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
};

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
  int nfaces(int &nel_surf,int &nel_vol) { nel_surf=4; nel_vol=8; return 24; }
  //@}
};

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Instantiation of the #Surf2Vol# class for
    quads with line surface elements. */
class Line2Quad : public Surf2Vol {
public:
  /// Constructor. Sets #use_exterior_normal_m#
  Quad2Hexa(const char *geom,int ndim,int nel,
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
  int nfaces(int &nel_surf,int &nel_vol) { nel_surf=2; nel_vol=4; return 4; }
  //@}
};
#endif

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Allows to compute integrals, or any kind of function
    that returns a set of scalar values over surface
    elements but depending on the gradient of state function
    normal to the surface. */
class embedded_gatherer : public Elemset { 
  /** Contains the interpolation/integration
      object (defines the element geometry). */
  Surf2Vol *sv_gp_data;
  /// Number of Gauss points. 
  int npg;
  int nel_surf, nel_vol, layer, layers;
public: 
  /** Constructor. Initializes the pointer 
      to interpolation/integration element. */
  embedded_gatherer() : sv_gp_data(NULL) { }
  /** Destructor. Safely destroys #sv_gp_data# because
      initialized to null pointer. */
  ~embedded_gatherer() { delete sv_gp_data; }
  /** Finds the volume element for each surface element and 
      rotates it in order to give the face to standard position. */
  void initialize();
  /** Number of elements in the global vector 
      that the element will contribute. */
  int gather_length;
  /** Loops over elements and Gauss points. Computes state vectors, 
      gradients, shape functions, and passes them to user as
      callback functions in order to compute contributions. */
  ASSEMBLE_FUNCTION;
  /// Only accepts the #gather# jobinfo. 
  ASK_FUNCTION;
  
  /** @name User defined callback functions
      @doc Deriving the class and defining these
      functions you can define several integrators.
  */
  //@{
  /** Number of surface nodes, desired number of
      nodes in the volume element. */
  virtual void surface_nodes(int &nel_surf,int &nel_vol)=0;
  /// Called \emph{before} the element loop. 
  virtual void init() { }

  /** Called for each element, \emph{before} the Gauss point loop. 
      @param k (input) element number in the elemeset. */ 
  virtual void element_hook(int k) {}

  /** Called at each Gauss point loop. Compute the contribution at
      each Gauss point. 
      @param pg_values (output) your contributions at the Gauss point
      @param u (input) the field state at time n+1
      @param uold (input) the field state at time n
      @param grad_u (input) the gradient of field state of variables at time n+1
      @param grad_uold (input) the gradient of field state of variables at time n
      @param xpg (input) coordinates of the Gauss point
      @param n (input) normal to the surface (if #ndimel=ndim-1#)
      between real and master coordinates
      @param wpgdet (input) the weight of the Gauss point
      @param time (input) time value ($t^\alpha$) */ 
  virtual void set_pg_values(vector<double> &pg_values,FastMat2 &u,
			     FastMat2 &uold,FastMat2 &grad_u, FastMat2 &grad_uold,
			     FastMat2 &xpg,FastMat2 &n, double wpgdet,double time)=0;

  /** Called \textbf{after} the element loop. May be used for
      clean-up operations. */
  virtual void clean() {};
  //@}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Computes the force (viscous and pressure) for Navier-Stokes eqs.
class visc_force_integrator : public embedded_gatherer { 
private:
  /// Auxiliary matrix, contains the force on the Gauss point. 
  FastMat2 force, moment, x_center, dx, strain_rate,
    sigma, Omega, rigid_grad_u;
  /// Flag to compute moments or not, number of dimensions
  int compute_moment, ndim_m;
  /// Viscosity
  double viscosity;
public:
  /// Get dimension number, define #compute_moment#, viscosity
  void init();
  /** Adds force (pressure+viscous) for a given Gauss point. 
      @param (input)
      @return a reference to the matrix.
  */ 
  void set_pg_values(vector<double> &pg_values,FastMat2 &u,
		     FastMat2 &uold,FastMat2 &grad_u, FastMat2 &grad_uold, 
		     FastMat2 &xpg,FastMat2 &n, double wpgdet,double time);
  /// Return parameters number of nodes on face, and on volume elements. 
  void surface_nodes(int &nel_surf,int &nel_vol);
};

#endif
