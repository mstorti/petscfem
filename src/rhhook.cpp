//__INSERT_LICENSE__
//$Id: rhhook.cpp,v 1.3 2006/04/11 16:44:46 mstorti Exp $

#include <src/debug.h>
#include <src/fem.h>
#include <src/readmesh.h>
#include <src/util2.h>
#include <src/util3.h>
#include <src/texthf.h>
#include <src/hook.h>
#include <src/rhhook.h>
#include <src/dvector.h>
#include <src/dvecpar.h>

extern int MY_RANK, SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_hfields_hook::
init(Mesh &mesh_a,Dofmap &dofmap_a,
     const char *name_a) {
  int ierr;

  mesh_p = &mesh_a;
  nu = mesh_a.nodedata->nu;
  ndim = mesh_a.nodedata->ndim;
  nnod = mesh_a.nodedata->nnod;

  //o File where to read the H fields
  TGETOPTDEF_S_ND(mesh_a.global_options,string,read_hfields_file,"Hfields.tmp");

  //o Ignore additional rows in the file (may be fictitious nodes)
  TGETOPTDEF_ND(mesh_a.global_options,int,read_hfields_ignore_extra_nodes,0);

  read_hfields();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_hfields_hook::
time_step_pre(double time,int step) {
  read_hfields();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_hfields_hook::read_hfields() {
  dvector<double> hfields;
  dvector_read_parallel(read_hfields_file.c_str(),hfields);

  int nrows, nH = nu-ndim;
  if (read_hfields_ignore_extra_nodes) {
    assert(hfields.size() >= nH*nnod);
    assert(hfields.size() % nH ==0);
    nrows = hfields.size()/nH;
  } else {
    assert(hfields.size() == nH*nnod);
    nrows = nnod;
  }
  hfields.reshape(2,nrows,nH);
  double *coords = mesh_p->nodedata->nodedata;
#define COORDSP(j,k) VEC2(coords,j,k,nu)
  for (int j=0; j<nnod; j++)
    for (int k=0; k<nH; k++) 
      COORDSP(j,ndim+k) = hfields.e(j,k);
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_hfields_hook::
time_step_post(double time,int step,
	       const vector<double> &gather_values) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void read_hfields_hook::close() {
}
#endif
