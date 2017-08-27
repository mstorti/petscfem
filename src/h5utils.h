// -*- mode: c++ -*-
//__INSERT_LICENSE__
#ifndef PETSCFEM_H5UTILS_H
#define PETSCFEM_H5UTILS_H

// HDF5 utils
void vec_gather(MPI_Comm comm,Vec v,vector<double> &values);
void h5petsc_mat_save(Mat J, const char *filename);
int h5petsc_vec_save(Vec x,const char *filename,const char *varname);
void h5_dvector_read(const char *file,const char *dset,dvector<double> &w);
void h5_dvector_write(dvector<double> &w,const char *filename,const char *varname);

#endif
