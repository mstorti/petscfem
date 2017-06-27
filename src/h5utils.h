// -*- mode: c++ -*-
//__INSERT_LICENSE__
#ifndef PETSCFEM_H5UTILS_H
#define PETSCFEM_H5UTILS_H

// HDF5 utils
void h5petsc_mat_save(Mat J, const char *filename);
int h5petsc_vec_save(Vec x,const char *filename,const char *varname);
void h5_dvector_read(const char *file,const char *dset,dvector<double> &w);

#endif
