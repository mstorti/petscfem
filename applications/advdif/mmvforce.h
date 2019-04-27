#ifndef PETSCFEM_MMVFORCE_H
#define PETSCFEM_MMVFORCE_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Forced mesh move hook for ALE.
// Moves the mesh with a function defined by the user. 
class mmv_force_hook_t : public Hook {
private:
  // File pattern to save the coordinates of the moving mesh
  string xalepattern;
  // Internal pointer to the coordinates of PF
  double *coords_buff;
  // Save coords to file frequency, current frame
  int nsaverot,frame;
  
public:
  // Reference coordinates. Coordinates at t{n} and t{n+1}
  dvector<double> coords0,coords;
  // Space dimension, number of nodes, number of columns in
  // coords array
  int ndim, nnod, nu;
  // Implemented hook functions
  mmv_force_hook_t();
  void init(Mesh &mesh_a,Dofmap &dofmap_a,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();

  // These are functions to be used by the user
  virtual void
  init_user(Mesh &mesh_a,Dofmap &dofmap_a,const char *name) { }
  virtual void
  time_step_post_user(double time,int step,
                      const vector<double> &gather_values) { }
  virtual void close_user() {}
  // Thi is the main function. It must compute the current
  // coordinates in columns [0,ndim) of COORDS. It can use
  // the coords at t{n} stored at columns [ndim,2*ndim) or
  // the initial coordinates in COORDS0
  virtual void
  compute_coords(dvector<double> &coords0,
                 dvector<double> &coords,
                 double time,double step) { }
};

#endif
