#ifndef PETSCFEM_MMVFORCE_H
#define PETSCFEM_MMVFORCE_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Forced mesh move hook for ALE.
// Moves the mesh with a function defined by the user. 
class mmv_force_hook_t : public Hook {
private:
  // File pattern to save the coordinates of the moving mesh
  string xalepattern;
  double *coords_buff;
  int nsaverot,frame;
  
public:
  dvector<double> coords0,coords;
  int ndim, nnod, nu;
  mmv_force_hook_t();
  void init(Mesh &mesh_a,Dofmap &dofmap_a,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();

  virtual void
  init_user(Mesh &mesh_a,Dofmap &dofmap_a,const char *name) { }
  virtual void
  time_step_pre_user(double time,int step) { } 
  virtual void
  time_step_post_user(double time,int step,
                      const vector<double> &gather_values) { }
  virtual void close_user() {}
  virtual void
  compute_coords(dvector<double> &coords0,
                 dvector<double> &coords,
                 double time,double step) { }
};

#endif
