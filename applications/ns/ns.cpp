//__INSERT_LICENSE__
//$Id: ns.cpp,v 1.20.2.1 2001/09/07 22:57:53 mstorti Exp $
 
#include <malloc.h>

#include "../../src/fem.h"
#include "../../src/readmesh.h"
#include "../../src/getprop.h"
#include "../../src/utils.h"
#include "../../src/util2.h"
#include "../../src/sttfilter.h"

#include "fracstep.h"
#include "nsi_tet.h"

static char help[] = "Basic finite element program.\n\n";

int MY_RANK,SIZE;
TextHashTable *GLOBAL_OPTIONS;

//debug:=
int TSTEP=0;

int print_internal_loop_conv_g=0;

     
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int MyKSPMonitor(KSP ksp,int n,double rnorm,void *dummy)
{
  int      ierr;

  if (print_internal_loop_conv_g) 
  PetscPrintf(PETSC_COMM_WORLD,
	      "iteration %d KSP Residual_norm = %14.12e \n",n,rnorm);
  return 0;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "bless_elemset"
void bless_elemset(char *type,Elemset *& elemset) {
  //  SET_ELEMSET_TYPE(internal)
  //    SET_ELEMSET_TYPE(fracstep) Por ahora lo desactivamos hasta que
  // hagamos la interfase
  // SET_ELEMSET_TYPE(nsi_tet)
    SET_ELEMSET_TYPE(nsi_tet_les)
    SET_ELEMSET_TYPE(nsi_tet_les_fm2)
    SET_ELEMSET_TYPE(nsi_tet_les_ther)
    SET_ELEMSET_TYPE(nsi_tet_keps)
    SET_ELEMSET_TYPE(bcconv_ns_fm2)
    SET_ELEMSET_TYPE(bcconv_nsther_fm2)
    SET_ELEMSET_TYPE(wall)
    SET_ELEMSET_TYPE(wallke)
	{
	printf("not known elemset \"type\": %s\n",type);
	exit(1);
	}
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** For a given state vector, prints the state for some nodes.
    @author M. Storti
    @param (input) filename file where to write the vector. May contain
    relative directories. 
    @param s (input) State vector
    @param dofmap (input) corresponding dofmap 
    @param node_list (input) set of nodes to print. 
*/ 
int print_some(const char *filename,const State &s,Dofmap *dofmap,
	       set<int> & node_list);
#endif


//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
#undef __FUNC__
#define __FUNC__ "main"
int main(int argc,char **args) {

  Vec     x, dx, xold,
    dx_step, res;		// approx solution, RHS, residual
  Mat     A_tet;		// linear system matrix 
  SLES    sles_tet;		// linear solver context
  PC      pc_tet;		// preconditioner context 
  KSP     ksp_tet;		// Krylov subspace method context
  double  norm, *sol, scal;	// norm of solution error
  int     ierr, i, n = 10, col[3], its, flg, size, node,
    jdof, k, kk, nfixa,
    kdof, ldof, lloc, nel, nen, neq, nu,
    myrank;
  // Initialize time
  Time time,time_star; time.set(0.);
  GlobParam glob_param;

  char fcase[FLEN+1];
  Dofmap *dofmap = new Dofmap;
  Mesh *mesh;
  arg_list argl;
  vector<double> hmin;
  hmin.resize(1);

  PetscInitialize(&argc,&args,(char *)0,help);
  print_copyright();

  // Start registering functions
  Amplitude::initialize_function_table();

  // Get MPI info
  MPI_Comm_size(PETSC_COMM_WORLD,&SIZE);
  MPI_Comm_rank(PETSC_COMM_WORLD,&MY_RANK);

  //  if (size != 1) SETERRA(1,0,"This is a uniprocessor example only!");
  ierr = OptionsGetString(PETSC_NULL,"-case",fcase,FLEN,&flg);
  CHKERRA(ierr);
  if (!flg) {
    PetscPrintf(PETSC_COMM_WORLD,
		"Option \"-case <filename>\" not passed to PETSc-FEM!!\n");
    PetscFinalize();
    exit(0);
  }

  // Read the mesh
  read_mesh(mesh,fcase,dofmap,neq,SIZE,MY_RANK);

  GLOBAL_OPTIONS = mesh->global_options;

  //o Dimension of the problem.
  GETOPTDEF(int,ndim,3);

  PetscPrintf(PETSC_COMM_WORLD,"After readmesh...\n");

  //o Number of inner iterations for the global non-linear
  // Newton  problem. 
  GETOPTDEF(int,nnwt,1);
  //o Tolerance to solve the non-linear system (global Newton).
  GETOPTDEF(double,tol_newton,1e-8);
  //o Update jacobian only until n-th Newton subiteration. 
  // Don't update if null. 
  GETOPTDEF(int,update_jacobian_iters,0);

  //o Absolute tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  GETOPTDEF(double,atol,1e-6);
  //o Relative tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  GETOPTDEF(double,rtol,1e-3);
  //o Divergence tolerance to solve the monolithic linear
  // system (Newton linear subiteration).
  GETOPTDEF(double,dtol,1e+3);
  //o Krylov space dimension in solving the monolithic linear
  // system (Newton linear subiteration) by GMRES.
  GETOPTDEF(int,Krylov_dim,50);
  //o Maximum iteration number in solving the monolithic linear
  // system (Newton linear subiteration).
  GETOPTDEF(int,maxits,Krylov_dim);
  //o Prints convergence in the solution of the GMRES iteration. 
  GETOPTDEF(int,print_internal_loop_conv,0);
  //o After computing the analytic Jacobian, Computes the
  // Jacobian in order to verify the analytic one. 
  GETOPTDEF(int,verify_jacobian_with_numerical_one,0);
  //o After computing the linear system solves it and prints Jacobian,
  // right hand side and solution vector, and stops. 
  GETOPTDEF(int,print_linear_system_and_stop,0);
  //o Measure performance of the 'comp\_mat\_res' jobinfo. 
  GETOPTDEF(int,measure_performance,0);

  print_internal_loop_conv_g=print_internal_loop_conv;

  //o Chooses the preconditioning operator. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,preco_type,jacobi);
  // I had to do this since `c_str()' returns `const char *'
  char *preco_type_ = new char[preco_type.size()+1];
  strcpy(preco_type_,preco_type.c_str());

  //o Sets the save frequency in iterations 
  GETOPTDEF(int,nsave,10);
  //o Sets the frequency save for the ``rotary save'' mechanism. 
  //i_tex ../../doc/nsdoc.tex rotary_save
  GETOPTDEF(int,nsaverot,100);
  //o Sets the number of states saved in a given file
  // in the ``rotary save'' mechanism (see \ref{sec:rotary_save}
  GETOPTDEF(int,nrec,1000000);
  //o Sets the number of files in the ``rotary save'' mechanism. 
  // (see \ref{sec:rotary_save})
  GETOPTDEF(int,nfile,1);

  //o Sets the save frequency in iterations for the ``print some''
  // mechanism. 
  //i_tex ../../doc/nsdoc.tex print_some
  GETOPTDEF(int,nsome,10000);
  //o The number of time steps. 
  GETOPTDEF(int,nstep,10000);
  //o The time step.
  GETOPTDEF(double,Dt,0.);
  //o Flag if steady solution or not (uses Dt=inf). If \verb+steady+
  // is set to 1, then the computations are as if $\Dt=\infty$. 
  // The value of \verb+Dt+ is used for printing etc... If \verb+Dt+
  // is not set and \verb+steady+ is set then \verb+Dt+ is set to one.
  GETOPTDEF(int,steady,0);
  if (steady && Dt==0.) Dt=1.;
  // Set values to be passed through global options
  glob_param.steady=steady;
  glob_param.Dt=Dt;
  //o Trapezoidal method parameter. \verb+alpha=1+:
  // Backward Euler. \verb+alpha=0+: Forward Euler.
  // \verb+alpha=0.5+: Crank-Nicholson. 
  GETOPTDEF(double,alpha,1.); 
  glob_param.alpha=alpha;
  //o Use the LES/Smagorinsky turbulence model. 
  GETOPTDEF(int,LES,0);

  //o The pattern to generate the file name to save in for
  // the rotary save mechanism.
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_pattern,outvector%d.out);

  //o The name of the file to save the state vector. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file,outvector.out);

  //o Name of file where to read the nodes for the ``print some'' 
  // feature. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,print_some_file,);
  //o Name of file where to save node values for the ``print some'' 
  // feature. 
  TGETOPTDEF_S(GLOBAL_OPTIONS,string,save_file_some,outvsome.out);
  //o Print, after execution, a report of the times a given option
  // was accessed. Useful for detecting if an option was used or not.
  GETOPTDEF(int,report_option_access,1);

  set<int> node_list;
  print_some_file_init(mesh->global_options,
		       print_some_file.c_str(),
		       save_file_some.c_str(),node_list);

  // initialize vectors
  dofmap->create_MPI_vector(x);
  State sx(x,time);		// convert to State

  LPFilterGroup filter(GLOBAL_OPTIONS,sx,Dt);

#if 0
  const int NT=200;
  const double DT=0.1, omega=1.;
  for (int kk=0; kk<NT; kk++) {
    double t=double(kk)*DT;
    double a=1.+sin(omega*t);
    sx.set_cnst(a);
    ff.update();
    PetscPrintf(PETSC_COMM_WORLD,"t = %f,  a = %f\n",t,a);
    //print_some("filter.some",ff.state(),dofmap,node_list);
    print_some("filter.some.tmp",sx,dofmap,node_list);
    print_some("filter.some.tmp",ff.state(),dofmap,node_list);
  }
  PetscFinalize();
  exit(0);
#endif
  
  ierr = VecDuplicate(x,&xold); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx_step); CHKERRA(ierr);
  ierr = VecDuplicate(x,&dx); CHKERRA(ierr);
  ierr = VecDuplicate(x,&res); CHKERRA(ierr);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // initialize state vectors
  scal=0;
  ierr = VecSet(&scal,x); CHKERRA(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"Before computing profile.\n");

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // COMPUTE ACTIVE PROFILE

  VOID_IT(argl);
  argl.arg_add(&A_tet,PROFILE);
  ierr = assemble(mesh,argl,dofmap,"comp_mat",&time); CHKERRA(ierr); 

#if 0 //dbg
  VOID_IT(argl);
  argl.arg_add(&A_tet,OUT_MATRIX);
  for (int jjj=0; jjj<10; jjj++) {
    printf("[loop iter %d]\n",jjj);
    ierr = assemble(mesh,argl,dofmap,"comp_mat",&time); CHKERRA(ierr); 
  }
  PetscFinalize();
  exit(0);
#endif

  // ierr =  MatSetOption(A_tet,MAT_NEW_NONZERO_ALLOCATION_ERR);
  PetscPrintf(PETSC_COMM_WORLD,"After computing profile.\n");

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Build octree for nearest neighbor calculation

  vector<double> *data_pts_ = new vector<double>;
  vector<ElemToPtr> *elemset_pointer = new vector<ElemToPtr>;
  WallData *wall_data;
  if (LES) {
    VOID_IT(argl);
    argl.arg_add(data_pts_,USER_DATA);
    argl.arg_add(elemset_pointer,USER_DATA);
    Elemset *elemset=NULL;
    argl.arg_add(elemset,USER_DATA);
    ierr = assemble(mesh,argl,dofmap,"build_nneighbor_tree",&time); CHKERRA(ierr); 
    PetscPrintf(PETSC_COMM_WORLD,"After nearest neighbor tree.\n");

    wall_data = new WallData(data_pts_,elemset_pointer,ndim);

    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // Find nearest neighbor for each volume element
    VOID_IT(argl);
    argl.arg_add(wall_data,USER_DATA);
    ierr = assemble(mesh,argl,dofmap,"get_nearest_wall_element",
		    &time); CHKERRA(ierr); 
  }

#if 0 
  ANNpoint point = annAllocPt(ndim);
  ANNidx nn_idx;
  ANNpoint nn = annAllocPt(ndim);
  ANNdist dist;
  int elem;


  while (1) {
    printf("\n\nEnter point: > ");
    scanf("%lf %lf %lf",&(point[0]),&(point[1]),&(point[2]));
    printf("\n");

    wall_data->nearest(point,elemset,elem,nn_idx,nn,dist);
    printf("Nearest neighbor id point %d, coords: %f %f %f, elem %d, elemset %p\n",
	   nn_idx,nn[0],nn[1],nn[2],elem,elemset);
  }
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  ierr = opt_read_vector(mesh,x,dofmap,MY_RANK); CHKERRA(ierr);

  // Filter *filter(x,*mesh);

  for (int tstep=1; tstep<=nstep; tstep++) {
    TSTEP=tstep; //debug:=
    time_star.set(time.time()+alpha*Dt);
    time.inc(Dt);
    if (print_internal_loop_conv_g) 
      PetscPrintf(PETSC_COMM_WORLD,
		  " --------------------------------------\n");
    PetscPrintf(PETSC_COMM_WORLD,"Time step: %d, time: %g\n",tstep,time.time());
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    // TET ALGORITHM

    // Inicializacion del paso
    ierr = VecCopy(x,dx_step);
    ierr = VecCopy(x,xold);

    double normres_external;
    for (int inwt=0; inwt<nnwt; inwt++) {

      // Initialize step
      int update_jacobian = !update_jacobian_iters || inwt<update_jacobian_iters;

      // Compute wall stresses
      VOID_IT(argl);
      argl.arg_add(&x,IN_VECTOR);
      argl.arg_add(&xold,IN_VECTOR);
      // wait_from_console("antes de comp_shear_vel"); 
      ierr = assemble(mesh,argl,dofmap,"comp_shear_vel",
		      &time_star); CHKERRA(ierr);
      // wait_from_console("despues de comp_shear_vel"); 

      // Communicate wall stresses amoung different processors
      // This is obsolete. Now we control this from the loop
      // over elements in the `wall' elemset. 

      // Mon Oct 23 19:13:39 ART 2000. Now I think that I need this
      // because the at each processor there is not a global version
      // of the state.
      if (LES && SIZE>1) {
	VOID_IT(argl);
	ierr = assemble(mesh,argl,dofmap,"communicate_shear_vel",
			&time_star); CHKERRA(ierr);
      }
      // wait_from_console("despues de communicate_shear_vel"); 

      scal=0;
      ierr = VecSet(&scal,res); CHKERRA(ierr);
      if (update_jacobian) {
	ierr = MatZeroEntries(A_tet); CHKERRA(ierr);
      }

      VOID_IT(argl);
      argl.arg_add(&x,IN_VECTOR);
      argl.arg_add(&xold,IN_VECTOR);
      argl.arg_add(&res,OUT_VECTOR);
      if (update_jacobian) argl.arg_add(&A_tet,OUT_MATRIX);
      argl.arg_add(&hmin,VECTOR_MIN);
      argl.arg_add(&glob_param,USER_DATA);
      argl.arg_add(wall_data,USER_DATA);

      const char *jobinfo = (update_jacobian ? "comp_mat_res" : "comp_res");

      // In order to measure performance
      if (measure_performance) {
	ierr = measure_performance_fun(mesh,argl,dofmap,jobinfo,
			&time_star); CHKERRA(ierr);
	PetscFinalize();
	exit(0);
      }

      ierr = assemble(mesh,argl,dofmap,jobinfo,&time_star); CHKERRA(ierr);

      Viewer matlab;
      if (verify_jacobian_with_numerical_one) {
	ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			       "output.m",&matlab); CHKERRA(ierr);
	ierr = ViewerSetFormat(matlab,
			       VIEWER_FORMAT_ASCII_MATLAB,"atet"); CHKERRA(ierr);
	ierr = MatView(A_tet,matlab);

	Mat A_tet_c;
	ierr = MatDuplicate(A_tet,MAT_DO_NOT_COPY_VALUES,&A_tet_c); CHKERRA(ierr);
	// ierr = MatCopy(A_tet,A_tet_c,SAME_NONZERO_PATTERN);

	ierr = MatZeroEntries(A_tet); CHKERRA(ierr);
	ierr = MatZeroEntries(A_tet_c); CHKERRA(ierr);

	VOID_IT(argl);
	argl.arg_add(&x,PERT_VECTOR);
	argl.arg_add(&xold,IN_VECTOR);
	argl.arg_add(&A_tet_c,OUT_MATRIX_FDJ);
	if (update_jacobian) argl.arg_add(&A_tet,OUT_MATRIX);
	argl.arg_add(&hmin,VECTOR_MIN);
	argl.arg_add(&glob_param,USER_DATA);
	argl.arg_add(wall_data,USER_DATA);
	ierr = assemble(mesh,argl,dofmap,jobinfo,
			&time_star); CHKERRA(ierr);

	ierr = ViewerSetFormat(matlab,
			       VIEWER_FORMAT_ASCII_MATLAB,"atet_fdj"); CHKERRA(ierr);
	ierr = MatView(A_tet_c,matlab);
	PetscFinalize();
	exit(0);
      }

      // SLES para el esquema global
//  	if (inwt>0) {
//  	  ierr = SLESDestroy(sles_tet); CHKERRA(ierr);  
//  	}
      ierr = SLESCreate(PETSC_COMM_WORLD,&sles_tet); CHKERRA(ierr);
      ierr = SLESSetOperators(sles_tet,A_tet,
			      A_tet,SAME_NONZERO_PATTERN); CHKERRA(ierr);
      ierr = SLESGetKSP(sles_tet,&ksp_tet); CHKERRA(ierr);
      ierr = SLESGetPC(sles_tet,&pc_tet); CHKERRA(ierr);

      ierr = KSPSetType(ksp_tet,KSPGMRES); CHKERRA(ierr);

	// NEW!!!!!!!! ====================================
      ierr = KSPSetPreconditionerSide(ksp_tet,PC_RIGHT);
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

	//ierr = KSPSetType(ksp_tet,KSPBICG); CHKERRA(ierr);
      ierr = KSPGMRESSetRestart(ksp_tet,Krylov_dim); CHKERRA(ierr);
      ierr = KSPSetTolerances(ksp_tet,rtol,atol,dtol,maxits);

      ierr = PCSetType(pc_tet,preco_type_); CHKERRA(ierr);
      ierr = KSPSetMonitor(ksp_tet,MyKSPMonitor,PETSC_NULL);

      ierr = SLESSolve(sles_tet,res,dx,&its); CHKERRA(ierr); 

      if (print_linear_system_and_stop) {
	ierr = ViewerASCIIOpen(PETSC_COMM_WORLD,
			       "system.dat",&matlab); CHKERRA(ierr);
	ierr = ViewerSetFormat(matlab,
			       VIEWER_FORMAT_ASCII_MATLAB,"atet"); CHKERRA(ierr);
	ierr =  SLESView(sles_tet,VIEWER_STDOUT_SELF);

	ierr = MatView(A_tet,matlab);

	ierr = ViewerSetFormat(matlab,
			       VIEWER_FORMAT_ASCII_MATLAB,"res"); CHKERRA(ierr);
	ierr = VecView(res,matlab);

	ierr = ViewerSetFormat(matlab,
			       VIEWER_FORMAT_ASCII_MATLAB,"dx"); CHKERRA(ierr);
	ierr = VecView(dx,matlab);
	
	PetscFinalize();
	exit(0);

      }	

      double normres;
      ierr  = VecNorm(res,NORM_2,&normres); CHKERRA(ierr);
      if (inwt==0) normres_external = normres;
      PetscPrintf(PETSC_COMM_WORLD,
		  "Newton subiter %d, norm_res  = %10.3e, update Jac. %d\n",
		  inwt,normres,update_jacobian);

      // update del subpaso
      scal= 1./alpha;
      ierr = VecAXPY(&scal,dx,x);

#if 0
      ierr = VecView(x,VIEWER_STDOUT_WORLD); CHKERRA(ierr);
      PetscFinalize();
      exit(0);
#endif

      ierr = SLESDestroy(sles_tet); CHKERRA(ierr);
    } // end of loop over Newton subiteration (inwt)

    // error difference
    scal = -1.0;
    ierr = VecAXPY(&scal,x,dx_step);
    ierr  = VecNorm(dx_step,NORM_2,&norm); CHKERRA(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"============= delta_u = %10.3e\n",norm);
    if (normres_external < tol_newton) {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Tolerance on newton loop reached:  || R ||_0,  norm_res =%g < tol = %g\n",
		  normres_external,tol_newton);
      break;
    }

    print_vector_rota(save_file_pattern.c_str(),x,dofmap,&time,
		      tstep-1,nsaverot,nrec,nfile);
  
    if (tstep % nsave == 0)
      print_vector(save_file.c_str(),x,dofmap,&time);

    filter.update(time);
    if (print_some_file!="" && tstep % nsome == 0) {
      print_some(save_file_some.c_str(),x,dofmap,node_list,&time);
      filter.print_some("filter.some.tmp",dofmap,node_list);
    }

  }

  if (report_option_access && MY_RANK==0) TextHashTable::print_stat();
  print_vector(save_file.c_str(),x,dofmap,&time);

  ierr = VecDestroy(x); CHKERRA(ierr); 
  ierr = VecDestroy(dx); CHKERRA(ierr); 
  ierr = VecDestroy(res); CHKERRA(ierr); 

  ierr = MatDestroy(A_tet); CHKERRA(ierr); 

#ifdef DEBUG_MALLOC_USE
  fclose(malloc_log);
#endif
  PetscFinalize();
  exit(0);
}

