//__INSERT_LICENSE__
//$Id: mmove.cpp,v 1.12 2002/12/02 21:51:44 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/texthf.h>
#ifdef USE_NEWMAT
#include <newmatap.h>
#endif

#include "nsi_tet.h"
#include "adaptor.h"
#include "mmove.h"

extern GlobParam *GLOB_PARAM;

void mesh_move::init() {

  int ierr;
  FastMat2 C;
  double c1 = sqrt(1./3.), c2 = sqrt(1./6.);
  // These are the gradient of shape functions for a master
  // tetra, with unit edge length with nodes at 
  // [+-1/2 0 0], [0,sqrt(3)/2,0], [0,1/sqrt(6),sqrt(2/3)]
  double c[12] = {-1., -c1, -c2, +1., -c1, -c2, 
		  0., 2*c1, -c2, 0., 0., 3*c2};
  
  G.resize(2,ndim,ndim); // the metric tensor
  xlocp.resize(1,nel*ndim); // The perturbed coordinates
  xloc0.resize(1,nel*ndim); // The perturbed coordinates
  assert(ndof==ndim);
  assert(ndim==2 || ndim==3);
  // assert(nel==ndim+1); // Only for triangles in 2D, tetras in 3D
  J.resize(2,ndim,ndim);
  dNdxi.resize(2,ndim,nel);

#ifdef USE_NEWMAT
  GG.ReSize(ndim);
  D.ReSize(ndim);
#endif

  if (ndim==2) {
    if (nel==3) {
      dNdxi.setel(-sin(M_PI/3)*cos(M_PI/6),1,1);
      dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,1);
      dNdxi.setel(+sin(M_PI/3)*cos(M_PI/6),1,2);
      dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,2);
      dNdxi.setel(0                       ,1,3);
      dNdxi.setel(+sin(M_PI/3)            ,2,3);
    } else if (nel==4) {
      double cq[] = {-1,-1,1,-1,1,1,-1,1};
      C.resize(2,nel,ndim).set(cq).t();
      dNdxi.set(C);
    } else PETSCFEM_ERROR("Only tringles ad quads in 2D: nel %d\n",nel);
  } else {
#if 0
    dNdxi.set(0.);
    for (int k=1; k<=3; k++) dNdxi.setel(1.,k,k);
    dNdxi.setel(-1.,1,4);
    dNdxi.setel(-1.,2,4);
    dNdxi.setel(-1.,3,4);
#endif
    
    C.resize(2,nel,ndim).set(c).t();
    dNdxi.set(C);
  }
  res_Dir.resize(2,nel,ndim);
  init_dfun();
  //o Perturbation scale length for increment in computing
  // the Jacobian with finite differences. 
  TGETOPTDEF_ND(thash,double,epsilon_x,1.e-4);
}

void mesh_move_eig::init_dfun() {
  int ierr;
  //o Distortion coefficient
  TGETOPTDEF_ND(thash,double,c_distor,1.);
  //o Functional exponent
  TGETOPTDEF_ND(thash,double,distor_exp,1.);
  //o Distortion coefficient
  TGETOPTDEF_ND(thash,double,c_volume,1.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double mesh_move::distor_fun(FastMat2 & xlocp) {
  xlocp.reshape(2,nel,ndim);
  J.prod(xlocp,dNdxi,-1,1,2,-1);
  xlocp.reshape(1,nel*ndim);
  G.prod(J,J,-1,1,-1,2);
  return distor_fun_G(G);
}

double mesh_move_eig::distor_fun_G(FastMat2 &G) {
  double df,la1,la2,la3,vol;
#ifdef USE_NEWMAT
  for (int i=1; i<=ndim; i++) {
    for (int j=1; j<=ndim; j++) {
      GG(i,j) = G.get(i,j);
    }
  }
  EigenValues(GG,D);
  la1 = D(1,1);
  la2 = D(2,2);
  if (ndim==3) la3 = D(3,3);
#else
  D.seig(G);
  la1=D.get(1);
  la2=D.get(2);
  if (ndim==3) la3=D.get(3);
#endif
#if 1
  double diffla;
  if (ndim==2) {
    vol = la1*la2;
    diffla = (la1-la2)*(la1-la2);
    diffla /= vol;
  } else if (ndim==3) {
    diffla = (la1-la2)*(la1-la2) + (la2-la3)*(la2-la3) + (la1-la3)*(la1-la3);
    vol = la1*la2*la3;
    diffla /= pow(vol,2./3.);
  } else assert(0);
  df = c_distor * pow(diffla,distor_exp) + c_volume * pow(vol,2.*distor_exp/double(ndim));
  return df;
#elif 0
  double p = 1.;
  double norm_D = D.norm_p_all(p);
  double norm_iD = D.norm_p_all(-p);
  return pow(norm_D * norm_iD, 1./p);
#else
  return D.max_all()/D.min_all();
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double mesh_move_rcond::distor_fun_G(FastMat2 & G) {
  iG.inv(G);
  double norm_G = G.max_abs_all();
  double norm_iG = iG.max_abs_all();
  return norm_G*norm_iG + 0.001*norm_G;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  double distor,distor_p,distor_m,
    distor_pp,distor_pm,distor_mp,distor_mm, d2f;
  double eps = epsilon_x;
  int nen = nel*ndim;
  xloc0.reshape(2,nel,ndim).set(xloc).add(state_new).reshape(1,nel*ndim);
  distor = distor_fun(xloc0);

  res.reshape(1,nel*ndim);
  for (int k=1; k<=nel*ndim; k++) {
    xlocp.set(xloc0);
    xlocp.setel(xloc0.get(k)+eps,k);
    distor_p = distor_fun(xlocp);

    xlocp.set(xloc0);
    xlocp.setel(xloc0.get(k)-eps,k);
    distor_m = distor_fun(xlocp);

    res.setel(-(distor_p-distor_m)/(2*eps),k);
  }
  res.reshape(2,nel,ndim);

  const FastMat2 *state = (glob_param->inwt==0 ? &state_old : &state_new);
  xloc0.reshape(2,nel,ndim).set(xloc).add(*state).reshape(1,nel*ndim);
  mat.reshape(2,nen,nen);
  for (int k=1; k<=nel*ndim; k++) {
    for (int l=k; l<=nel*ndim; l++) {
      xlocp.set(xloc0);
      xlocp.setel(xloc0.get(k)+eps,k);
      xlocp.setel(xlocp.get(l)+eps,l);
      distor_pp = distor_fun(xlocp);
      
      xlocp.set(xloc0);
      xlocp.setel(xloc0.get(k)+eps,k);
      xlocp.setel(xlocp.get(l)-eps,l);
      distor_pm = distor_fun(xlocp);
      
      xlocp.set(xloc0);
      xlocp.setel(xloc0.get(k)-eps,k);
      xlocp.setel(xlocp.get(l)+eps,l);
      distor_mp = distor_fun(xlocp);
      
      xlocp.set(xloc0);
      xlocp.setel(xloc0.get(k)-eps,k);
      xlocp.setel(xlocp.get(l)-eps,l);
      distor_mm = distor_fun(xlocp);

      d2f = (distor_pp - distor_pm - distor_mp + distor_mm)/(4.*eps*eps);
      mat.setel(d2f,k,l);
    }
  }

  for (int k=2; k<=nel*ndim; k++) {
    for (int l=1; l<=k-1; l++) {
      d2f = mat.get(l,k);
      mat.setel(d2f,k,l);
    }
  }

  mat.reshape(4,nel,ndim,nel,ndim);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int mmove_hook::write_mesh(const State &s,const char *filename,
			   const int append=0) {
#define XNOD(j,k) VEC2(xnod,j,k,dofmap->ndof)
  int myrank;
  double *vseq_vals,*sstate,*xnod;
  Vec vseq;
  const Vec &x = s.v();

  MPI_Comm_rank(PETSC_COMM_WORLD,&myrank);

  // fixme:= Now we can make this without a scatter. We can use
  // the version of get_nodal_value() with ghost_values. 
  int neql = (myrank==0 ? dofmap->neq : 0);
  int ierr = VecCreateSeq(PETSC_COMM_SELF,neql,&vseq);  CHKERRQ(ierr);
  ierr = VecScatterBegin(x,vseq,INSERT_VALUES,
			 SCATTER_FORWARD,*dofmap->scatter_print); CHKERRA(ierr); 
  ierr = VecScatterEnd(x,vseq,INSERT_VALUES,
		       SCATTER_FORWARD,*dofmap->scatter_print); CHKERRA(ierr); 
  ierr = VecGetArray(vseq,&vseq_vals); CHKERRQ(ierr);

  xnod = mesh->nodedata->nodedata;

  if (myrank==0) {
    printf("Writing vector to file \"%s\"\n",filename);
    FILE *output;
    output = fopen(filename,(append == 0 ? "w" : "a" ) );
    if (output==NULL) {
      printf("Couldn't open output file\n");
      // fixme:= esto esta mal. Todos los procesadores
      // tienen que llamar a PetscFinalize()
      exit(1);
    }

    int ndof=dofmap->ndof;
    double dval;
    const Time *t = &s.t();
    for (int k=1; k<=dofmap->nnod; k++) {
      for (int kldof=1; kldof<=ndof; kldof++) {
	dofmap->get_nodal_value(k,kldof,vseq_vals,t,dval);
	fprintf(output,"%12.10e  ",XNOD(k-1,kldof-1)+dval);
      }
      fprintf(output,"\n");
    }
    fclose(output);
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mmove_hook::time_step_post(double time,int step,
				const vector<double> &gather_values) {
  int ierr = write_mesh(*GLOB_PARAM->state,"remeshing.dat",1); 
  assert(ierr=0);
}

