//__INSERT_LICENSE__
//$Id merge-with-petsc-233-55-g52bd457 Fri Oct 26 13:57:07 2007 -0300$
 
#include <src/fem.h>
#include <src/utils.h>
#include <src/utils.h>

#ifdef USE_HDF5
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
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int h5petsc_vec_save(Vec x,const char *filename,const char *varname) {
  PetscObjectSetName((PetscObject)x,varname);
  PetscViewer viewer;
  int ierr;
  ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,filename,
			     FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);
  ierr = VecView(x,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  return ierr;
}
#else
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void h5petsc_mat_save(Mat J, const char *filename) {
  PETSCFEM_ERROR0("Not compiled with HDF5 support");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int h5petsc_vec_save(Vec x,const char *filename,const char *varname) {
  PETSCFEM_ERROR0("Not compiled with HDF5 support");
}
#endif
