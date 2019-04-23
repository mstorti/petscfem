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
// #include <applications/advdif/mmvforce.h>

using namespace std;

static char help[] = "Basic finite element program.\n\n";

extern int print_internal_loop_conv_g,
  consistent_supg_matrix_g,
  local_time_step_g,
  comp_mat_each_time_step_g,
  verify_jacobian_with_numerical_one;

#define VECVIEW(name,label) \
ierr = PetscViewerSetFormat(matlab, \
		       PETSC_VIEWER_ASCII_MATLAB,#label); \
ierr = VecView(name,matlab); CHKERRA(ierr)

// PETSc now doesn't have the string argument that
// represents the variable name so that I will use this
// wrapper until I find how to set names in Ascii matlab
// viewers.
#define PetscViewerSetFormat_WRAPPER(viewer,format,name) \
          PetscViewerSetFormat(viewer,format)

Hook *advdif_hook_factory(const char *name);
extern Mesh *GLOBAL_MESH;
Dofmap *GLOBAL_DOFMAP;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
bool ajk_comp(ajk_t a,ajk_t b) {
  if (a.j!=b.j) return a.j<b.j;
  return a.k<b.k;
}

// The Chimera hook that specializes the code for a
// particular case
chimera_hook_t *CHIMERA_HOOK_P=NULL;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
int chimera_mat_shell_t::init0(Mat A_) {

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
int chimera_mat_shell_t::init1(Vec x,Vec res) {

  // List of nodes at the boundaries of W1 and W2 (includes
  // external and internal boundaries).
  // Call a hook from the user
  CHIMERA_HOOK_P->mark_bdry_nodes(ebdry,ibdry);

  // Replace all the rows for the external and internal
  // boundary nodes for the identity matrix Replace the rows
  set<int> fixed = ebdry;
  for (auto &q : ibdry) fixed.insert(q);
  rows.clear();
  for (auto &node : fixed) {
    int jeq = node2eq[node];
    rows.push_back(jeq);
  }
  fixed.clear();
  printf("Imposed rows (external+internal) bdries %zu\n",rows.size());
  
  // Load the interpolators (computed in Octave right now probably)
  dvector<double> w;
  h5_dvector_read("./interp.h5","/z/value",w);
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

  int ierr;
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
int mat_mult(Mat Ashell,Vec x,Vec y) {
  void *ctx;
  int ierr = MatShellGetContext(Ashell,&ctx); CHKERRQ(ierr);
  chimera_mat_shell_t &cms = *(chimera_mat_shell_t*)ctx;
  // Make the base contribution to the matrix-vector product
  ierr = MatMult(cms.A,x,y); CHKERRQ(ierr);
  // Add the extra therm (restrictions between fos by chimera)
  cms.mat_mult(x,y);
  return 0;
}

// FIXME:= is it needed??
void init_hooks();

//-------<*>-------<*>-------<*>-------<*>-------<*>-------
#undef __FUNC__
#define __FUNC__ "chimera_main"
int chimera_main() {
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
  // This is just because otherwise the main program doesn't
  // link the functions. FIXME:= is it needed??
  init_hooks();
  
  PetscBool flg;
  int ierr;
  Vec     x, dx, xold, res; /* approx solution, RHS, residual*/
  PFMat *A,*AA;			// linear system matrix
  PFMat *A_tet, *A_tet_c;
  double  *sol, scal, normres, normres_ext=NAN;    /* norm of solution error */
  int     i, n = 10, col[3], its, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, ndim, nel, nen, neq, nu,
    myrank;
  // nu:= dimension of the state vector per node
  PetscScalar  neg_one = -1.0, one = 1.0, value[3];
  PetscScalar *px;
  char fcase[FLEN+1];
  Darray *da; // este me parece que se puede sacar!!
  //Elemset *elemset;
  Dofmap *dofmap = NULL;
  dofmap = new Dofmap;
  Mesh *mesh = NULL;
  vector<double> dtmin(1,0.);
  Vec a;
  GlobParam glob_param;
  GLOB_PARAM = &glob_param;
  string save_file_res;

  print_copyright();
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "-------- Generic Advective-Diffusive  module ---------\n"
	      "-------- CHIMERA MODULE ---------\n");

  Debug debug(0,PETSCFEM_COMM_WORLD);
  GLOBAL_DEBUG = &debug;

  ierr = PetscOptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg); CHKERRA(ierr);
  if (!flg) {
    PetscPrintf(PETSCFEM_COMM_WORLD,
		"Option \"-case <filename>\" not passed to PETSc-FEM!!\n");
    PetscFinalize();
    exit(0);
  }

  // Read data
  ierr = read_mesh(mesh,fcase,dofmap,neq,SIZE,MY_RANK); CHKERRA(ierr);
  GLOBAL_OPTIONS = mesh->global_options;
  GLOBAL_DOFMAP = dofmap;

  //o Activate debugging
  GETOPTDEF(int,activate_debug,0);
  if (activate_debug) {
    debug.activate();
    Debug::init();
  }
  //o Activate printing in debugging
  GETOPTDEF(int,activate_debug_print,0);
  if (activate_debug_print) debug.activate("print");
  //o Activate report of memory usage
  GETOPTDEF(int,activate_debug_memory_usage,0);
  if (activate_debug_memory_usage) debug.activate("memory_usage");

  //o Absolute tolerance when solving a consistent matrix
  GETOPTDEF(double,atol,1e-6);
  //o Relative tolerance when solving a consistent matrix
  GETOPTDEF(double,rtol,1e-3);
  //o Divergence tolerance when solving a consistent matrix
  GETOPTDEF(double,dtol,1e+3);
  //o Maximum iterations when solving a consistent matrix
  GETOPTDEF(int,maxits,150);
  //o Prints the convergence history when solving a consistent matrix
  GETOPTDEF(int,print_internal_loop_conv,0);
  print_internal_loop_conv_g=print_internal_loop_conv;
  //o Measure performance of the  #comp_mat_res#  jobinfo.
  GETOPTDEF(int,measure_performance,0);

  //o Save state vector frequency (in steps)
  GETOPTDEF(int,nsave,10);
  //o Save state vector frequency with the ``rotary save''
  // mechanism. (see \ref{sec:rotary_save})
  GETOPTDEF(int,nsaverot,100);
  //o Sets the number of states saved in a given file
  // in the ``rotary save'' mechanism (see \ref{sec:rotary_save}
  GETOPTDEF(int,nrec,1000000);
  //o Sets the number of files in the ``rotary save'' mechanism.
  // (see \ref{sec:rotary_save})
  GETOPTDEF(int,nfile,1);

  //o The number of time steps.
  GETOPTDEF(int,nstep,10000);
  //o Output CPU time statistics for frequency in time steps.
  GETOPTDEF(int,nstep_cpu_stat,10);
  //o After computing the linear system prints Jacobian and
  // right hand side and stops..
  GETOPTDEF(int,print_linear_system_and_stop,0);
  //o Solve system before  #print\_linear_system_and_stop#
  GETOPTDEF(int,solve_system,1);
  //o If #print_linear_system_and_stop# is active,
  // then print system in this Newton iteration
  GETOPTDEF(int,inwt_stop,0);
  //o If #print_linear_system_and_stop# is active,
  // then print system in this time step
  GETOPTDEF(int,time_step_stop,1);
  //o Print the residual each  #nsave#  steps.
  GETOPTDEF(int,print_residual,0);

  //o Sets the save frequency in iterations for the ``print some''
  // mechanism. (see doc in the Navier-Stokes module)
  GETOPTDEF(int,nsome,10000);
  //o Name of file where to read the nodes for the ``print some''
  // feature.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,print_some_file,);
  //o Name of file where to save node values for the ``print some''
  // feature.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_some,outvsome.out);
  //o Access mode to the ``some'' file. If 0 rewind file. If 1
  //  append to previous  results.
  TGETOPTDEF(GLOBAL_OPTIONS,int,save_file_some_append,1);
  //o Print, after execution, a report of the times a given option
  // was accessed. Useful for detecting if an option was used or not.
  GETOPTDEF(int,report_option_access,1);
  //o Update jacobian each $n$-th time step.
  GETOPTDEF(int,update_jacobian_steps,0);
  //o Use IISD (Interface Iterative Subdomain Direct) or not.
  GETOPTDEF(int,use_iisd,0);
  //o Type of solver. May be  #iisd#  or  #petsc# .
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,solver,petsc);
  if (use_iisd) solver = string("iisd");

  // Use IISD (Interface Iterative Subdomain Direct) or not.
  // A_tet = (use_iisd ? &IISD_A_tet : &PETSc_A_tet);
  A  = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());
  AA = PFMat::dispatch(dofmap->neq,*dofmap,solver.c_str());

  set<int> node_list;
  print_some_file_init(mesh->global_options,
		       print_some_file.c_str(),
		       save_file_some.c_str(),node_list,
		       save_file_some_append);

  // warning: passed to advective.cpp via a global variable
  //o Uses consistent SUPG matrix for the temporal term or not.
  GETOPTDEF(int,consistent_supg_matrix,0);
  consistent_supg_matrix_g = consistent_supg_matrix;
  //o Chooses automatically the time step from the
  // selected Courant number
  GETOPTDEF(int,auto_time_step,1);
  // warning: passed to advective.cpp via a global variable
  //o Chooses a time step that varies locally. (Only makes sense
  // when looking for steady state solutions.
  GETOPTDEF(int,local_time_step,1);
  local_time_step_g=local_time_step;
  //o The Courant number.
  GETOPTDEF(double,Courant,0.6);
  //o Time step.
  GETOPTDEF(double,Dt,0.);
  //o Flag if steady solution or not (uses Dt=inf). If  #steady#
  // is set to 1, then the computations are as if $\Dt=\infty$.
  // The value of  #Dt#  is used for printing etc... If  #Dt#
  // is not set and  #steady#  is set then  #Dt#  is set to one.
  GETOPTDEF(int,steady,0);
  if (steady && Dt==0.) Dt=1.;
  glob_param.Dt = Dt;
  glob_param.steady = steady;
  //o The parameter of the trapezoidal rule
  // for temporal integration.
  GETOPTDEF(double,alpha,1.);
  glob_param.alpha=alpha;
#define ALPHA (glob_param.alpha)
  if (ALPHA>0.) {
    consistent_supg_matrix = 1;
    auto_time_step=0;		// Don't know how to do auto_time_step
				// in the implicit case...
    local_time_step=0;
  }
  //o Number of iterations in the Newton loop. (
  // for the implicit method.
  GETOPTDEF(int,nnwt,3);
  if (ALPHA==0.) nnwt=1;

  //o Flag for launching RENORM process
  GETOPTDEF(int,RENORM_flag,0);

  comp_mat_each_time_step_g =
    consistent_supg_matrix || local_time_step;

  //o Counts time from here.
  GETOPTDEF(double,start_comp_time,0.);
  //o Counts time from here. (superseded by
  //  #start_comp_time# for compatibility with Navier-Stokes
  //  module).
  GETOPTDEF(double,start_time,NAN);
  if (!ISNAN(start_time) && start_comp_time!=0.0)
    start_comp_time = start_time;

  //o Tolerance when solving with the mass matrix.
  GETOPTDEF(double,tol_mass,1e-3);
  //o Tolerance when solving the sublinear problem
  // at each iteration.
  GETOPTDEF(double,tol_linear,1e-3);
  //o Tolerance when solving the non-linear problem
  // for the implicit case.
  GETOPTDEF(double,tol_newton,1e-3);
  //o Tolerance when solving for a steady state
  GETOPTDEF(double,tol_steady,0.);
  //o Relaxation factor for the Newton iteration
  GETOPTDEF(double,omega_newton,1.);
  //o Computes jacobian of residuals and prints to a file.
  //  May serve to debug computation of the analytic jacobians.
  TGETOPTDEF_ND(mesh->global_options,int,verify_jacobian_with_numerical_one,0);
  //o Check whether the states are finite (not Inf or NaN).
  //  If this happens the program stops.
  GETOPTDEF(int,check_for_inf,0);

#define INF INT_MAX
  //o Update jacobian each $n$-th time step.
  GETOPTDEF(int,update_jacobian_start_steps,INF);
  //o Update jacobian only until n-th Newton subiteration.
  // Don't update if null.
  GETOPTDEF(int,update_jacobian_iters,1);
  PETSCFEM_ASSERT0(update_jacobian_iters>=1,"Out of range");
  //o Update jacobian each $n$-th Newton iteration
  GETOPTDEF(int,update_jacobian_start_iters,INF);
  PETSCFEM_ASSERT0(update_jacobian_start_iters>=0,"Out of range");
#undef INF

  vector<double> gather_values;
  //o Number of ``gathered'' quantities.
  GETOPTDEF(int,ngather,0);
  //o Print values in this file
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,gather_file,gather.out);
  // Initialize gather_file
  FILE *gather_file_f;
  if (MY_RANK==0 && ngather>0) {
    gather_file_f = fopen(gather_file.c_str(),"w");
    //fprintf(gather_file_f,"");
    fclose(gather_file_f);
  }

  //o Chooses the preconditioning operator.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,preco_type,jacobi);
  // I had to do this since `c_str()' returns `const char *'
  char *preco_type_ = new char[preco_type.size()+1];
  strcpy(preco_type_,preco_type.c_str());

  Time time,time_star;
  time.set(start_comp_time);
  glob_param.time = &time;

  //o The pattern to generate the file name to save in for
  // the rotary save mechanism.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_pattern,outvector%d.out);

  //o Filename for saving the state vector.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file,outvector.out);
  save_file_res = save_file + string(".res");

  //o Use HDF5 for saving linear system
  TGETOPTDEF(mesh->global_options,int,use_hdf5,0);

#if 0
  PetscViewer matlab;
  ierr = PetscViewerASCIIOpen(PETSCFEM_COMM_WORLD,
			 "matns.m",&matlab); CHKERRA(ierr);
#endif
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:

  dofmap->create_MPI_vector(x);

  // initialize vectors
  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);

  // Set pointers in glob_param
  glob_param.x = x;
  glob_param.xold = xold;
  State state(x,time),state_old(xold,time);
  glob_param.state = &state;
  glob_param.state_old = &state_old;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // initialize state vectors
  scal=0;
  ierr = VecSet(x,scal); CHKERRA(ierr);

  arg_list argl,arglf;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Hook stuff
  HookList hook_list;
  hook_list.init(*mesh,*dofmap,advdif_hook_factory);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Compute  profiles
  debug.trace("Computing profiles...");
  VOID_IT(argl);
  argl.arg_add(A,PROFILE|PFMAT);
  ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr);
  debug.trace("After computing profile.");

  if (ADVDIF_CHECK_JAC) {
    VOID_IT(argl);
    argl.arg_add(AA,PROFILE|PFMAT);
    ierr = assemble(mesh,argl,dofmap,"comp_prof",&time); CHKERRA(ierr);
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  ierr = opt_read_vector(mesh,x,dofmap,MY_RANK); CHKERRA(ierr);

  // Initialize the CHIMERA stuff
  Mat Ashell;
  chimera_mat_shell_t cms;

  Mat Apetsc = A->get_petsc_mat();
  MatSetOption(Apetsc,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
  cms.init0(Apetsc);
  int nlocal;
  // ierr = VecGetSize(dx,&neq);CHKERRQ(ierr);
  ierr = VecGetLocalSize(dx,&nlocal);CHKERRQ(ierr);
  ierr = MatCreateShell(PETSC_COMM_WORLD,nlocal,nlocal,neq,neq,
                        &cms,&Ashell);
  MatShellSetOperation(Ashell,MATOP_MULT,(void (*)(void))(&mat_mult));
  KSP ksp;         /* linear solver context */
  PC pc;           /* preconditioner context */
  ierr = KSPCreate(PETSCFEM_COMM_WORLD,&ksp);CHKERRQ(ierr);

  ierr = KSPSetOperators(ksp,Ashell,cms.A,
                         DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,
                          PETSC_DEFAULT); CHKERRQ(ierr);
  ierr = KSPMonitorSet(ksp,KSPMonitorDefault,NULL,NULL); CHKERRQ(ierr);
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // This is for taking statistics of the
  // CPU time consumed by a time steptime
  Chrono chrono;
#define STAT_STEPS 5
  double cpu[STAT_STEPS],cpuav;
  int update_jacobian_this_step,update_jacobian_this_iter;
  for (int tstep=1; tstep<=nstep; tstep++) {
    time_star.set(time.time()+alpha*Dt);
    // Take cputime statistics
    if (MY_RANK==0) {
      if (tstep>1) {
	int last=(tstep-2) % STAT_STEPS;
	cpu[last] = chrono.elapsed();
	if (tstep==2) {
	  for (int jstep=1; jstep<STAT_STEPS; jstep++) {
	    cpu[jstep] = cpu[0];
	  }
	}
	cpuav = 0;
	for (int jstep=0; jstep<STAT_STEPS; jstep++)
	  cpuav+= cpu[jstep];
	cpuav /= STAT_STEPS;
	if (tstep % nstep_cpu_stat == 0)
	  printf("CPU per iter: %f, average [over last %d]: %f\n",
		 cpu[last],STAT_STEPS,cpuav);
      }
    }
    chrono.start();

    if (RENORM_flag)
      ierr = read_vector("state-ren.tmp",x,dofmap,MY_RANK); CHKERRA(ierr);

    ierr = VecCopy(x,xold);
    //    hook_list.time_step_pre(time_star.time(),tstep);
    hook_list.time_step_pre(time.time()+Dt,tstep); //hook needs t_{n+1}

    for (int inwt=0; inwt<nnwt; inwt++) {
      // Initializes res
      scal=0;
      ierr = VecSet(res,scal); CHKERRA(ierr);

      if (comp_mat_each_time_step_g) {

	//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
	// ierr = A->build_ksp(GLOBAL_OPTIONS); CHKERRA(ierr);

	ierr = A->clean_mat(); CHKERRA(ierr);
        if (ADVDIF_CHECK_JAC) {
          ierr = AA->clean_mat(); CHKERRA(ierr);
        }
	VOID_IT(argl);
	argl.arg_add(&xold,IN_VECTOR);

        if (!ADVDIF_CHECK_JAC) argl.arg_add(&x,IN_VECTOR);
        else argl.arg_add(&x,PERT_VECTOR);

	argl.arg_add(&res,OUT_VECTOR);
	argl.arg_add(&dtmin,VECTOR_MIN);
	argl.arg_add(A,OUT_MATRIX|PFMAT);
	argl.arg_add(&glob_param,USER_DATA);
        if (ADVDIF_CHECK_JAC)
          argl.arg_add(AA,OUT_MATRIX_FDJ|PFMAT);

	if (measure_performance) {
	  ierr = measure_performance_fun(mesh,argl,dofmap,"comp_res",
					 &time_star); CHKERRA(ierr);
	  PetscFinalize();
	  exit(0);
	}
	debug.trace("Before residual computation...");
	ierr = assemble(mesh,argl,dofmap,"comp_res",&time_star); CHKERRA(ierr);
	debug.trace("After residual computation.");

	if (!print_linear_system_and_stop || solve_system) {
	  debug.trace("Before solving linear system...");
          cms.init1(x,res);
          ierr = KSPSolve(ksp,res,dx); CHKERRQ(ierr);
#if 0          
          DBG_MM = 1;
          ierr = MatMult(Ashell,dx,res); CHKERRQ(ierr);
          DBG_MM = 0;
#endif
	  debug.trace("After solving linear system.");
	}
      } else {

	VOID_IT(argl);
	argl.arg_add(&x,IN_VECTOR);
	argl.arg_add(&res,OUT_VECTOR);
	argl.arg_add(&dtmin,VECTOR_MIN);

	if (measure_performance) {
	  ierr = measure_performance_fun(mesh,argl,dofmap,"comp_res",
					 &time_star); CHKERRA(ierr);
	  PetscFinalize();
	  exit(0);
	}
	ierr = assemble(mesh,argl,dofmap,"comp_res",&time_star); CHKERRA(ierr);

	if (!print_linear_system_and_stop || solve_system) {
	  ierr = A->solve(res,dx); CHKERRA(ierr);
	}
      }

      if (print_linear_system_and_stop &&
	  inwt==inwt_stop && tstep==time_step_stop) {
        if (!use_hdf5) {
          PetscPrintf(PETSCFEM_COMM_WORLD,
                      "Printing residual and matrix for debugging and stopping..\n");
          PetscViewer matlab;
          ierr = PetscViewerASCIIOpen(PETSCFEM_COMM_WORLD,
                                      "mat.output",&matlab); CHKERRA(ierr);
          ierr = PetscViewerSetFormat_WRAPPER(matlab,
                                              PETSC_VIEWER_ASCII_MATLAB,"res");

          PetscObjectSetName((PetscObject)res,"Vec_0");
          ierr = VecView(res,matlab);
          if (solve_system) {
            ierr = PetscViewerSetFormat_WRAPPER(matlab,
                                                PETSC_VIEWER_ASCII_MATLAB,"dx");
            ierr = VecView(dx,matlab);
          }

          ierr = PetscViewerSetFormat_WRAPPER(matlab,
                                              PETSC_VIEWER_ASCII_MATLAB,"A");
          Mat AP = A->get_petsc_mat();
          PetscObjectSetName((PetscObject)AP,"Mat_1");
          ierr = A->view(matlab);
          print_vector(save_file_res.c_str(),res,dofmap,&time); // debug:=
          if (ADVDIF_CHECK_JAC) {
            ierr = PetscViewerSetFormat_WRAPPER(matlab,
                                                PETSC_VIEWER_ASCII_MATLAB,"AA");
            ierr = AA->view(matlab);
          }
        } else {
#ifdef USE_HDF5
          Mat AA = A->get_petsc_mat();
          h5petsc_mat_save(AA,"advdifsys.h5");
          h5petsc_vec_save(res,"res.h5","res");
#else
          PETSCFEM_ERROR0("Save in HDF5 format requested but code "
                         "was not compiled with HDF5 support\n");

#endif
        }
	PetscFinalize();
	exit(0);
      }

      ierr  = VecNorm(res,NORM_2,&normres); CHKERRA(ierr);
      if (check_for_inf) {
#if 0
        PETSCFEM_ASSERT0(VecIsFinite(res),
                         "Detected Inf or NaN values in residual vector");
#else
        PETSCFEM_ASSERT0(isfinite(normres),
                         "Detected Inf or NaN values in residual vector");
#endif
      }
      if (inwt==0) normres_ext = normres;
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  "Newton subiter %d, norm_res  = %10.3e\n",
		  inwt,normres);
      scal=omega_newton/alpha;
      ierr = VecAXPY(x,scal,dx);
      if (check_for_inf) {
#if 0
        PETSCFEM_ASSERT0(VecIsFinite(x),
                         "Detected Inf or NaN values in state vector");
#else
        double normx;
        ierr  = VecNorm(x,NORM_2,&normx); CHKERRA(ierr);
        PETSCFEM_ASSERT0(isfinite(normx),
                         "Detected Inf or NaN values in state vector");
#endif
      }
      if (normres < tol_newton) break;
      } // end of Newton loop

    // Prints residual and mass matrix in Matlab format
    // Define time step depending on strategy. Automatic time step,
    // local time step, etc...
    if (auto_time_step) Dt = Courant*dtmin[0];
    if (local_time_step) Dt = Courant;
    if (Dt<=0.) {
      PFEMERRQ("error: Dt<=0. You have to set Dt>0 or use "
	       "the auto_time_step flag.\n");
    }
    // SHV(Dt);
    time.inc(Dt);

    // Upgrade state vector.

    VOID_IT(argl);
    argl.arg_add(&x,IN_OUT_VECTOR);
    argl.arg_add(&xold,IN_VECTOR);
    ierr = assemble(mesh,argl,dofmap,"absorb_bc_proj",&time_star); CHKERRA(ierr);

    double time_=time;

    // Reuse dx for computing the delta state vector after re-projecting
    double delta_u;
    ierr = VecCopy(x,dx);
    scal=-1.;
    ierr = VecAXPY(dx,scal,xold);
    ierr  = VecNorm(dx,NORM_2,&delta_u); CHKERRA(ierr);

    if (tstep % nsave == 0){
      PetscPrintf(PETSCFEM_COMM_WORLD,
		  " --------------------------------------\n"
		  "Time step: %d\n"
		  " --------------------------------------\n",
		  tstep);

      print_vector(save_file.c_str(),x,dofmap,&time);
      if (print_residual)
	print_vector(save_file_res.c_str(),res,dofmap,&time);
    }

    if (ngather>0) {
      gather_values.resize(ngather,0.);
      for (int j=0; j<ngather; j++) gather_values[j] = 0.;
      arglf.clear();
      arglf.arg_add(&x,IN_VECTOR);
      arglf.arg_add(&xold,IN_VECTOR);
      arglf.arg_add(&gather_values,VECTOR_ADD);
      ierr = assemble(mesh,arglf,dofmap,"gather",&time);
      CHKERRA(ierr);
    }

    //    hook_list.time_step_post(time_star.time(),tstep,gather_values);
    hook_list.time_step_post(time.time()+Dt,tstep,gather_values);

    if (ngather>0) {
      // Print gathered values
      if (MY_RANK==0) {
	if (gather_file == "") {
	  printf("Gather results: \n");
	  for (unsigned int j=0; j < gather_values.size(); j++)
	    printf("v_component_%d = %12.10e\n",j,gather_values[j]);
	} else {
	  gather_file_f = fopen(gather_file.c_str(),"a");
	  for (unsigned int j=0; j<gather_values.size(); j++)
	    fprintf(gather_file_f,"%12.10e ",gather_values[j]);
	  fprintf(gather_file_f,"\n");
	  fclose(gather_file_f);
	}
      }
    }

    PetscPrintf(PETSCFEM_COMM_WORLD,
		"time_step %d, time: %g, delta_u = %10.3e\n",
		tstep,time_,delta_u);

    print_vector_rota(save_file_pattern.c_str(),x,dofmap,
		      &time,tstep-1,nsaverot,nrec,nfile);

    if (print_some_file!="" && tstep % nsome == 0)
      print_some(save_file_some.c_str(),x,dofmap,node_list,&time);

    if (normres_ext < tol_steady) break;

  }
  hook_list.close();

  print_vector(save_file.c_str(),x,dofmap,&time);
  if (print_residual)
    print_vector(save_file_res.c_str(),res,dofmap,&time);
  if (report_option_access && MY_RANK==0) TextHashTable::print_stat();

  ierr = VecDestroy(&x); CHKERRA(ierr);
  ierr = VecDestroy(&xold); CHKERRA(ierr);
  ierr = VecDestroy(&dx); CHKERRA(ierr);
  ierr = VecDestroy(&res); CHKERRA(ierr);
#ifdef DIAG_MAT_MATRIX
  ierr = MatDestroy(&A); CHKERRA(ierr);
#endif

  // ierr = VecZeroEntries(res);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&Ashell);CHKERRQ(ierr);

  delete A;
  delete AA;
  DELETE_SCLR(dofmap);
  DELETE_SCLR(mesh);

#ifdef DO_SIZE_STATS
  prod2_subcache_t::report_stats();
#endif

  PetscFinalize();
  exit(0);
}
