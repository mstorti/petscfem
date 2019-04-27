#include <src/debug.h>
#include <time.h>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <json/json.h>

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/pfmat.h>
#include <src/hook.h>
#include <src/h5utils.h>
#include <src/dvector.h>
#include <applications/advdif/chimera.h>
#include <applications/advdif/advective.h>
#include <ANN/ANN.h>
#include <tools/project/project.h>

using namespace std;
extern Mesh *GLOBAL_MESH;
extern Dofmap *GLOBAL_DOFMAP;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
bool ajk_comp(ajk_t a,ajk_t b) {
  if (a.j!=b.j) return a.j<b.j;
  return a.k<b.k;
}

// The Chimera hook that specializes the code for a
// particular case
chimera_hook_t *CHIMERA_HOOK_P=NULL;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int chimera_mat_shell_t::init(Mat A_) {

#ifdef USE_JSONCPP
  // Read data needed from a JSON file
  ifstream in("data.json");
  in >> opts;
  cout << "Input opts: ====================" << endl
       << opts << endl;
  nnod1 = opts["nnod1"].asInt();
  nnod2 = opts["nnod2"].asInt();
  nelem1 = opts["nelem1"].asInt();
  nelem2 = opts["nelem2"].asInt();
  printf("nnod1 %d, nelem1 %d, nnod2 %d, nelem2 %d\n",
         nnod1,nelem1,nnod2,nelem2);

  // Get the coordinates of the nodes
  xnod = GLOBAL_MESH->nodedata->nodedata;
  nu = GLOBAL_MESH->nodedata->nu;
  nnod = GLOBAL_MESH->nodedata->nnod;

#else
  PETSCFEM_ERROR0("Not compiled with JSONCPP support\n");
#endif
  
  // We prepare the system to solve A\res
  int ierr;
  // Store a pointer to the underlying PETSc matrix
  A = A_;
  // Check that only one processor is being used
  int size;
  MPI_Comm_size(PETSCFEM_COMM_WORLD,&size);
  PETSCFEM_ASSERT0(size==1,"Only one processor so far");

  // Get the dofmap in order to map equations to nodes
  Dofmap *dofmap = GLOBAL_DOFMAP;
  int nnod = dofmap->nnod;

  int neq = dofmap->neq;
  // So far only used for scalar problems (ndof==1)
  PETSCFEM_ASSERT0(dofmap->ndof==1,"Only 1 dof/node so far");
  // The RHS vector. If we will replace the equation for
  // some node JEQ to PHI[JEQ]=VAL then we have to set the
  // RHS[JEQ] to VAL and set the corresponding row to 1
  // (Identity) In this stage we just set the rows of the
  // matrix to Identity and the RHS.
  // COUNT will be the number of rows that are set
  int count=0;
  PETSCFEM_ASSERT0(nnod==neq,
                   "Not allowed Dirichlet conditions "
                   "through FIXA for PF-CHIMERA");
  for (int node=0; node<nnod; node++) {
    int m;
    const int *dofs;
    const double *coefs;
    // Dofmap works with base 1 nodes and dofs...
    dofmap->get_row(node+1,1,m,&dofs,&coefs);
    PETSCFEM_ASSERT0(m==1,"Dofmap is not a permutation of the identity!!");
    double tol=1e-5;
    PETSCFEM_ASSERT0(coefs[0]==1.0,"Dofmap is not identity!!");
    int jeq = dofs[0]-1;
    eq2node[jeq] = node;
    node2eq[node] = jeq;
  }
  return 0;
}

int DBG_MM=0;

#define XNOD(j,k) VEC2(xnod,j,k,nu)
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int chimera_mat_shell_t
::before_solve(Vec x,Vec res,double time,int step) {

  int ierr;
  // List of nodes at the boundaries of W1 and W2 (includes
  // external and internal boundaries).
  // Call a hook from the user
  CHIMERA_HOOK_P->mark_bdry_nodes(ebdry,ibdry,time,step);

  // Replace all the rows for the external and internal
  // boundary nodes for the identity matrix Replace the rows
  set<int> fixed = ebdry;
  for (auto &q : ibdry) fixed.insert(q);
  // This stores the rows that are fixed
  vector<int> rows;
  for (auto &node : fixed) {
    int jeq = node2eq[node];
    rows.push_back(jeq);
  }
  fixed.clear();
  printf("Imposed rows (external+internal) bdries %zu\n",rows.size());

  TGETOPTDEF(GLOBAL_OPTIONS,int,use_octave_interpolator,0);
  Nodedata *nd = GLOBAL_MESH->nodedata;
  double *xnod = nd->nodedata;
  int nnod = nd->nnod;
  int nu = nd->nu;
  TGETOPTDEF(GLOBAL_OPTIONS,int,ndim,0);
  assert(ndim>0);

  PETSCFEM_ASSERT0(nu==2*ndim,"Not correct dims");

  // Get the PF current coords
  dvector<double> xale;
  xale.a_resize(2,nnod,ndim);
  for (int j=0; j<nnod; j++)
    for (int k=0; k<ndim; k++)
      xale.e(j,k) = XNOD(j,k);
  
  dvector<double> w;
  if (use_octave_interpolator) {
    // Write coordinates to a HDF5 file so as to call then
    // an Octave scrpt and compute externally the
    // interpolators
    // Save the file
    const char *fname = "./coords.h5";
    if (!access(fname,F_OK)) unlink(fname);
    h5_dvector_write(xale,fname,"xale");
    // FIXME:= the computation of the interpolators should
    // go to the Chimera module
    system("octave-cli -qH mkint.m > mkinterpolators.log");
    // Load the interpolators (computed in Octave right now probably)
    h5_dvector_read("./interp.h5","/z/value",w);
  } else {
    // We extract the data (xnod,icone) from the PF internal
    // data.
    // We assume that we have only one elemset
    PETSCFEM_ASSERT0(Elemset::elemset_table.size()==1,
                     "More than one elemset defined. Not supported yet.");
    
    map<string,Elemset *>::iterator 
      q = Elemset::elemset_table.begin();
    Elemset *elemset = q->second;
    int nelem = elemset->size();
    int nel,ndof,nelprops;
    elemset->elem_params(nel,ndof,nelprops);
#define ICONE(j,k) VEC2(elemset->icone,j,k,nel)
    dvector<int> icone;
    icone.a_resize(2,nelem,nel);
    for (int e=0; e<nelem; e++) 
      for (int k=0; k<nel; k++)
        icone.e(e,k) = ICONE(e,k);
    printf("read %d elems nel=%d from PF Elemset %s\n",
           nelem,nel,elemset->name());

    // Number of neighbors to be search by ANN
    int nngbr=10;
    FemInterp fem_interp;
    fem_interp.print_area_coords="USE_RETVAL";
    fem_interp.init(10,ndof,ndim,xale,icone);
    dvector<double> u1;
    u1.a_resize(2,nnod,ndof);
    double z=0.0;
    u1.set(z);
    fem_interp.interp(xale,u1,w);
    printf("interpolator size: %d %d\n",w.size(0),w.size(1));
    exit(0);
  }
  xale.clear();

  int ncoef = w.size(0);
  printf("Loaded interpolators. ncoef %d\n",ncoef);
  PETSCFEM_ASSERT0(w.size(1)==3,"Bad z column size");
  
  z.clear();
  for (int l=0; l<ncoef; l++) {
    int
      j=dbl2int(w.e(l,0)),
      k=dbl2int(w.e(l,1));
    double a = w.e(l,2);
    ajk_t ajk(j,k,a);
    z.push_back(ajk);
  }
  w.clear();

  sort(z.begin(),z.end(),ajk_comp);
  zptr.clear();
  zptr.resize(nnod+1,-1);
  int jlast=0;
  for (int l=0; l<ncoef; l++) {
    int j = z[l].j,k = z[l].k;
    while (jlast<=j) {
      zptr[jlast++] = l;
    }
  }
  while (jlast<=nnod) zptr[jlast++] = ncoef;

  double *resp,*xp;
  ierr = VecGetArray(res,&resp); CHKERRQ(ierr);
  ierr = VecGetArray(x,&xp); CHKERRQ(ierr);
  // For external bdry nodes the RHS of the eq is simply 0,
  // because we assume homogeneous Dirichlet condition
  for (auto &jeq : ebdry) resp[jeq] = 0.0;
  // For internal bdry nodes we must set the difference
  // between the value of PHI and the interpolated value
  // from the other domain. But as we solve in incremental
  // form PHI = X+DX, so we have to put in the RHS
  // rhs{j} = phi{j} - sum{k} a{jk}*phi{k}
  // rhs{j} = x{j} - sum{k} a{jk}*x{k} + dx{j} - sum{k} a{jk}*dx{k}
  // The term dx{j}-sum{k} a{jk}*dx{k} added by
  // the MatMult product. 
  // So we put in the RHS the following:
  // rhs{j} = x{j} - sum{k} a{jk}*x{k}
  for (auto &node1 : ibdry) {
    resp[node1] = -xp[node1];
    int
      rstart = zptr[node1],
      rend = zptr[node1+1];
    PETSCFEM_ASSERT(rend>rstart,
                    "Can't find interpolator for boundary "
                    "node %d, x(%f,%f)",node1,XNOD(node1,0),XNOD(node1,1));
    double val=0.0,sumcoef=0.0;
    for (int l=rstart; l<rend; l++) {
      ajk_t &a = z[l];
      int node2=a.k;
      int jeqk = node2eq[a.k];
      if (DBG_MM && node1==0)
        printf("node2 %d, x %f %f, phi %f\n",
               node2,XNOD(node2,0),XNOD(node2,1),xp[jeqk]);
      sumcoef += a.ajk;
      val += a.ajk*xp[jeqk];
    }
    resp[node1] += val;
  }

  // Set the rows for the external and internal boundaries
  // to the identity
  int nrows = rows.size();
  ierr = MatZeroRows(A,nrows,rows.data(),1.0,NULL,NULL); CHKERRQ(ierr);
  ierr = VecRestoreArray(res,&resp); CHKERRQ(ierr);
  ierr = VecRestoreArray(x,&xp); CHKERRQ(ierr);
  
  return 0;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int chimera_mat_shell_t::mat_mult(Vec x,Vec y) {
  // Here we can add some extra contribution to the
  // residuals
  int ierr;
  double *xp,*yp;
  ierr = VecGetArray(x,&xp); CHKERRQ(ierr);
  ierr = VecGetArray(y,&yp); CHKERRQ(ierr);

  int j2 = 297;
  if (DBG_MM) printf("x[0] %f, x[%d] %f\n",xp[0],j2,xp[j2]);
  
  // Interpolate de values at the internal W1 bdry
  for (auto &node1 : ibdry) {
    int
      rstart = zptr[node1],
      rend = zptr[node1+1];
    PETSCFEM_ASSERT(rend>rstart,
                    "Can't find interpolator for boundary "
                    "node %d, x(%f,%f)",node1,XNOD(node1,0),XNOD(node1,1));
    // printf("node1 %d x(%f %f)\n",
    //        node1,XNOD(node1,0),XNOD(node1,1));
    // Look for node in the interpolator with largest coefficient
    double val=0.0,sumcoef=0.0;
    for (int l=rstart; l<rend; l++) {
      ajk_t &a = z[l];
      int node2=a.k;
      // printf("-> node2 %d x(%f %f) coef %g\n",
      //        node2,XNOD(node2,0),XNOD(node2,1),a.ajk);
      int jeqk = node2eq[a.k];
      if (DBG_MM && node1==0)
        printf("node2 %d, x %f %f, phi %f\n",
               node2,XNOD(node2,0),XNOD(node2,1),xp[jeqk]);
      sumcoef += a.ajk;
      val += a.ajk*xp[jeqk];
    }
    int jeq1 = node2eq[node1];
    if (DBG_MM && node1==0)
      printf("node1 %d x(%f %f) id %f, interp %f,sumcoef %f\n",
             node1,XNOD(node1,0),XNOD(node1,1),yp[jeq1],val,sumcoef);
    yp[jeq1] += -val;
  }
  if (DBG_MM) printf("y[0] %f\n",yp[0]);
  ierr = VecRestoreArray(x,&xp); CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&yp); CHKERRQ(ierr);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int chimera_mat_mult(Mat Ashell,Vec x,Vec y) {
  void *ctx;
  int ierr = MatShellGetContext(Ashell,&ctx); CHKERRQ(ierr);
  chimera_mat_shell_t &cms = *(chimera_mat_shell_t*)ctx;
  // Make the base contribution to the matrix-vector product
  ierr = MatMult(cms.A,x,y); CHKERRQ(ierr);
  // Add the extra therm (restrictions between fos by chimera)
  cms.mat_mult(x,y);
  return 0;
}
