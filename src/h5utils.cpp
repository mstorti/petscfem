#include <src/fem.h>
#include <src/utils.h>
#include <src/h5utils.h>

// Given a distributed PETSc vector `vec' gather all the
// ranges in all processor in a full vector of doubles
// `values'.  This full vector is available in all the
// processors.
// WARNING: this may be inefficient and non
// scalable, it is just a dirty trick to have access to all
// the values of the vectors in all the processor.
// Usage:
// Vec v;
// // ... create and fill v eith values at each processor
// // ... do the Assembly
// vector<double> values;
// vec_gather(MPI_COMM_WORLD,v,values);
// //... now you have all the elements of `v' in `values'
void vec_gather(MPI_Comm comm,Vec v,vector<double> &values) {
  // n: global size of vector
  // nlocal: local (PETSc) size
  int n,nlocal;
  // Get the global size
  VecGetSize(v,&n);
  // Resize the local buffer
  values.clear();
  values.resize(n,0.0);
  // Get the local size
  VecGetLocalSize(v,&nlocal);

  // Gather all the local sizes in order to compute the
  // counts and displs for the Allgatherv
  int size, myrank;
  MPI_Comm_rank(comm,&myrank);
  MPI_Comm_size(comm,&size);
  vector<int> counts(size),displs(size);
  MPI_Allgather(&nlocal,1,MPI_INT,
                &counts[0],1,MPI_INT,comm);
  displs[0]=0;
  for (int j=1; j<size; j++)
    displs[j] = displs[j-1] + counts[j-1];

  // Get the internal values of the PETSc vector
  double *vp;
  VecGetArray(v,&vp);
  // Do the Allgatherv to the local vector
  MPI_Allgatherv(vp,nlocal,MPI_DOUBLE,
                 &values[0],&counts[0],&displs[0],MPI_DOUBLE,comm);
  // Restore the array
  VecRestoreArray(v,&vp);
}

#ifdef USE_HDF5
#include "H5Cpp.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5petsc_mat_save(Mat J, const char *filename) {
  int m,n;
  MatGetSize(J,&m,&n);

  vector<int> rows,cols;
  vector<double> vals;
  int ncols;
  const int *colsp;
  const double *valsp;
  for (int j=0; j<m; j++) {
    MatGetRow(J,j,&ncols,&colsp,&valsp);
    for (int l=0; l<ncols; l++) {
      rows.push_back(j);
      cols.push_back(colsp[l]);
      vals.push_back(valsp[l]);
    }
  }

  H5::H5File file(filename,H5F_ACC_TRUNC);

  hsize_t ncoefs = rows.size();
  H5::DataSpace dataspace(1,&ncoefs);

  // Create the dataset.
  H5::DataSet dsrows =
    file.createDataSet("rows",H5::PredType::NATIVE_INT,dataspace);
  dsrows.write(rows.data(),H5::PredType::NATIVE_INT);

  H5::DataSet dscols =
    file.createDataSet("cols",H5::PredType::NATIVE_INT,dataspace);
  dscols.write(cols.data(),H5::PredType::NATIVE_INT);

  H5::DataSet dsvals =
    file.createDataSet("vals",H5::PredType::NATIVE_DOUBLE,dataspace);
  dsvals.write(vals.data(),H5::PredType::NATIVE_DOUBLE);

  vector<int> sz(2);
  sz[0] = m;
  sz[1] = n;
  hsize_t ssz = 2;
  H5::DataSpace ds2(1,&ssz);
  H5::DataSet dssz =
    file.createDataSet("sizes",H5::PredType::NATIVE_INT,ds2);
  dssz.write(sz.data(),H5::PredType::NATIVE_INT);
  file.close();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int h5petsc_vec_save(Vec x,const char *filename,
                     const char *varname) {
  vector<double> vx;
  vec_gather(PETSC_COMM_WORLD,x,vx);
  H5::H5File file(filename,H5F_ACC_TRUNC);
  hsize_t n = vx.size();
  H5::DataSpace dataspace(1,&n);
  // Create the dataset.
  H5::DataSet xdset =
    file.createDataSet(varname,H5::PredType::NATIVE_DOUBLE,dataspace);
  xdset.write(vx.data(),H5::PredType::NATIVE_DOUBLE);
  file.close();
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Split the file/dataset combination FDNAME in the FILE and
// DSET parts. For instance if fdname=./mesh.h5:/u/value then
// file=./mesh.h5, dset=/u/value
static void h5_file_dset_split(string &fdname,string &file,
                               string &dset) {
  int len = fdname.size();
  int j;
  for (j=0; j<len; j++) if (fdname[j]==':') break;
  PETSCFEM_ASSERT(j<len,"Couldn't find colon in HDF5 "
                  "file/dataset name: %s",fdname.c_str());  
  file = fdname.substr(0,j);
  dset = fdname.substr(j+1,len);
  // printf("fdname %s, file %s, dset %s\n",
  //        fdname.c_str(),file.c_str(),dset.c_str());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
template<class T>
void h5_dvector_read(const char *filename,
                     const char *dsetname,
                     dvector<T> &w,H5::PredType type) {
  if (!filename) return;
  H5::H5File h5file(filename,H5F_ACC_RDONLY);
  H5::DataSet dataset = h5file.openDataSet(dsetname);
  H5::DataSpace space = dataset.getSpace();
  int rank = space.getSimpleExtentNdims();
  assert(rank==2 || rank==1);

  vector<hsize_t> dims(rank);
  int ndims = space.getSimpleExtentDims(dims.data(),NULL);
  int sz = 1;
  for (int j=0; j<ndims; j++) sz *= dims[j];
  if (w.size()==0) w.resize(sz);
  PETSCFEM_ASSERT(sz==w.size(),"Wrong size, sz(w) %d, sz(h5) %d",
                  w.size(),sz);
  w.defrag();
  dataset.read(w.buff(),type);
  vector<int> shape(rank);
  for (int j=0; j<rank; j++) shape[j] = int(dims[j]);
  w.reshape(shape);
  h5file.close();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
template<class T>
void h5_dvector_read2(const char *fdname,
                     dvector<T> &w) {
  string file,dset, fds=fdname;
  h5_file_dset_split(fds,file,dset);
  h5_dvector_read(file.c_str(),dset.c_str(),w);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_read(const char *fdname,
                     dvector<double> &w) {
  h5_dvector_read2(fdname,w);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_read(const char *fdname,
                     dvector<int> &w) {
  h5_dvector_read2(fdname,w);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_read(const char *filename,
                     const char *dsetname,
                     dvector<double> &w) {
  h5_dvector_read(filename,dsetname,w,
                  H5::PredType::NATIVE_DOUBLE);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_read(const char *filename,
                     const char *dsetname,
                     dvector<int> &w) {
  h5_dvector_read(filename,dsetname,w,
                  H5::PredType::NATIVE_INT);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_write(dvector<double> &w,const char *filename,
                     const char *varname) {
  w.defrag();
  H5::H5File *filep=NULL;
  if (!access(filename,F_OK)) 
    filep = new H5::H5File(filename,H5F_ACC_RDWR);
  else 
    filep = new H5::H5File(filename,H5F_ACC_TRUNC);
  vector<int> shape;
  w.get_shape(shape);
  int rank=shape.size();
  PETSCFEM_ASSERT(rank>=0 && rank<=2,
                  "Not implemented yet rank %d",rank);  
  vector<hsize_t> hshape;
  if (rank>0) {
    for (int j=0; j<rank; j++) hshape.push_back(shape[j]);
  } else {
    rank=1;
    hshape.push_back(w.size());
  }
  H5::DataSpace dataspace(rank,hshape.data());
  // Create the dataset.
  H5::DataSet xdset =
    filep->createDataSet(varname,H5::PredType::NATIVE_DOUBLE,dataspace);
  xdset.write(w.buff(),H5::PredType::NATIVE_DOUBLE);
  filep->close();
  delete filep;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Converts double to int, checking that the double is really an int
int dbl2int(double z) {
  int k = int(z);
  PETSCFEM_ASSERT(fabs(z-double(k))==0.0,
                  "Double is not integer! %g",z);
  return k;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*> Read a
// Read a dvector stored as a double in HDF5 and convert to
// INT (checking that the values are really integer)
void h5_dvector_read_d2i(const char *file,const char *dset,
                         dvector<int> &w) {
  dvector<double> z;
  h5_dvector_read(file,dset,z);
  vector<int> shape;
  z.get_shape(shape);
  int n = z.size();
  w.resize(n);
  for (int j=0; j<n; j++)
    w.ref(j) = dbl2int(z.ref(j));
  z.reshape(shape);
}

#else
#define H5ERR PETSCFEM_ERROR0("Not compiled with HDF5 support")
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5petsc_mat_save(Mat J, const char *filename) { H5ERR; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int h5petsc_vec_save(Vec x,const char *filename,
                     const char *varname) { H5ERR; return 0; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_read(const char *filename,const char *dsetname,
                     dvector<double> &w) { H5ERR; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int h5petsc_vec_save(Vec x,const char *filename,
                     const char *varname) { H5ERR; return 0; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
template<class T>
void h5_dvector_read(const char *filename,
                     const char *dsetname,
                     dvector<T> &w,H5::PredType type) { H5ERR; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
template<class T>
void h5_dvector_read2(const char *fdname,
                     dvector<T> &w) { H5ERR; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_read(const char *fdname,
                     dvector<double> &w) { H5ERR; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_read(const char *fdname,
                     dvector<int> &w) { H5ERR; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_read(const char *filename,
                     const char *dsetname,
                     dvector<double> &w) { H5ERR; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_read(const char *filename,
                     const char *dsetname,
                     dvector<int> &w) { H5ERR; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5_dvector_write(dvector<double> &w,const char *filename,
                     const char *varname) { H5ERR; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*> Read a
void h5_dvector_read_d2i(const char *file,const char *dset,
                         dvector<int> &w) { H5ERR; }

#endif

