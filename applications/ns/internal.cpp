/*__INSERT_LICENSE__*/
//$Id: internal.cpp,v 1.5 2001/05/30 03:58:44 mstorti Exp $
 
#include "fem.h"
#include "readmesh.h"
#include "fastmat.h"
#include "getprop.h"
#include "internal.h"

#define MAXPROP 100


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "internal::assemble"
int internal::assemble(double *retval,Nodedata *nodedata,double *locst,
		       Dofmap *dofmap,int ijob,int myrank,
		       int el_start,int el_last,int iter_mode) {

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define RETVAL(iele,j,k) VEC3(retval,iele,j,nel,k,ndof)
#define RETVALMAT(iele,j,k,p,q) VEC5(retval,iele,j,nel,k,ndof,p,nel,q,ndof)

  int ierr=0;

  if (ijob==COMP_MAT) {
    printf("This elemset routine no longer computes matrices\n");
    exit(1);
  }

  if (!(ijob==COMP_VEC || ijob==COMP_FDJ)) {
    printf("Don't know how to compute ijob: %d\n",ijob);
    exit(1);
  }
    
  //#define NODEDATA(j,k) (nodedata[ndim*(j)+(k)])
  //#define NODEDATA(j,k) VEC2(nodedata,j,k,ndim)
#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 
#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define LOAD(j,k) (load[ndof*(j)+(k)]) 
#define IDENT(j,k) (ident[ndof*(j)+(k)]) 
#define JDOFLOC(j,k) VEC2(jdofloc,j,k,ndof)
  
  int locdof,kldof,lldof;
  char *value;

  // Unpack Elemset
  int npg,ndim;
  ierr = get_int(thash,"npg",&npg); CHKERRA(ierr);
  ierr = get_int(thash,"ndim",&ndim); CHKERRA(ierr);
  int nen = nel*ndof;

  // Unpack Dofmap
  int *ident,neq,nnod;
  //  double *load;
  ident = dofmap->ident;
  neq = dofmap->neq;
  nnod = dofmap->nnod;
  //  load = dofmap->load;
  // fixme:= esto por ahora lo dejo asi. Que eventualmente pueda
  // ser diferente el ndof de 
  //  ndof = dofmap->ndof;

  // Unpack nodedata
  int nu=nodedata->nu;
  if(nnod!=nodedata->nnod) {
    printf("nnod from dofmap and nodedata don't coincide\n");
    exit(1);
  }

  // allocate local vecs
  int kdof;
  //jdofloc = new int[nel*ndof];
//      Matrix matloc(nen,nen), xloc(nel,ndim), veccontr(nel,ndof),
//        locstate(nel,ndof);
  FastMat matloc(nen,nen), xloc(nel,ndim), veccontr(nel,ndof),
    locstate(nel,ndof),tmp(nel,ndof);

  nen = nel*ndof;

  // Physical properties
  int iprop=0, elprpsindx[MAXPROP]; double propel[MAXPROP];

#if 0
  void *dummy;                
  g_hash_table_foreach (elem_prop_names,&print_prop_hash_entry,dummy);
  exit(0);
#endif 

DEFPROP(conductivity)
#define COND (propel+conductivity_indx)
 
DEFPROP(propa)
#define PROPA (propel+propa_indx)
 
DEFPROP(propb)
#define PROPB (propel+propb_indx)
 
DEFPROP(propc)
#define PROPC (propel+propc_indx)
 
DEFPROP(propp)
#define PROPP (propel+propp_indx)

DEFPROP(propq)
#define PROPQ (propel+propq_indx)

  int nprops=iprop;
 
  Matrix vel(ndof,ndim); vel=0.;
  ierr = get_double(thash,"velocity",vel.Store(),1,ndof*ndim); CHKERRA(ierr);

  double kond;
  ierr = get_double(thash,"conductivity",&kond); CHKERRA(ierr);

  // Gauss Point data
  //  Matrix *shape, *dshapexi, *dshapex;
  //  double *wpg;
  char *geom;
  thash->get_entry("geometry",geom);
  GPdata gp_data(geom,ndim,nel,npg);

  // Definiciones para descargar el lazo interno
//    Matrix grad_fi, Jaco, iJaco, Uintri;
//    RowVector Pert, W;
//    //  LogAndSign L;
  double detJaco, UU, u2, Peclet, psi, tau;
  int elem, ipg,node, jdim, kloc,lloc,ldof;

  // make fastmat copy version
  FastMat **FM_shape,**FM_dshapexi;
  FM_shape = new (FastMat *)[npg];
  FM_dshapexi = new (FastMat *)[npg];
  for (ipg=0; ipg<npg; ipg++) {
    FM_shape   [ipg] = new FastMat(gp_data.shape[ipg]);
    FM_dshapexi[ipg] = new FastMat(gp_data.dshapexi[ipg]);
  }
  FastMat dshapex(ndim,nel),Jaco(ndim,ndim),iJaco(ndim,ndim),
    grad_fi(ndim,ndof),dshapext(nel,ndim);
  
  for (int k=el_start; k<=el_last; k++) {
    if (epart[k] != myrank+1) continue;
    load_props(propel,elprpsindx,nprops,&(ELEMPROPS(k,0)));
    //    printf("element %d, prop %f\n",k,ELEMPROPS(k,0));
    // Load local node coordinates in local vector
    for (kloc=0; kloc<nel; kloc++) {
      node = ICONE(k,kloc);
      for (jdim=0; jdim<ndim; jdim++) {
	xloc.set(kloc,jdim,NODEDATA(node-1,jdim));
      }
    }
    locstate.set(&(LOCST(k,0,0)));
    
    veccontr.set(0.);
#define SHAPE    (*FM_shape[ipg])
#define DSHAPEXI (*FM_dshapexi[ipg])
#define WPG      (gp_data.wpg[ipg])

    // loop over Gauss points
    for (ipg=0; ipg<npg; ipg++) {

      //      Matrix xpg = SHAPE * xloc;

      ierr = FMp(Jaco,DSHAPEXI,xloc); FMCHK;

      elem = k;

      ierr = FMdet(&detJaco,Jaco); FMCHK;
      
      if (detJaco <= 0.) {
	cout << "Jacobian of element " << elem << " is negative or null\n"
	     << " Jacobian: " << detJaco << endl ;
	exit(1);
      }
      ierr = FMinv(iJaco,Jaco); FMCHK;
      ierr = FMp(dshapex,iJaco,DSHAPEXI); FMCHK;

      /*-------<*>-------<*>-------<*>-------<*>-------<*>-------
	SUPG stuff
	---------<*>-------<*>-------<*>-------<*>-------<*>-----*/ 
      // SUPG perturbation function

      // Uintri := velocity vector in master element covariant coordinates 
      // u2 := square of the norm of the velocity vector 
      // Peclet := non-dimensional Peclet number
      // psi := magic function psi(Peclet)
      // tau := intrinsic time for the element
      // Pert := SUPG perturbation function (\tilde P)
      // W := SUPG weight function (N + \tilde P)
      // UUU:= is an estimation of u/h
      // fixme: this is not completely invariant for triangles,
      //        since the element master is not equilateral. 

      ierr = FMp(grad_fi,dshapex,locstate); FMCHK;
      ierr = FMtr(dshapext,dshapex); FMCHK;

//        for (int kdof=0; kdof<ndof; kdof++) {
//   	Uintri = iJaco * (vel.Row(kdof+1)).t();
//   	UU = sqrt(Uintri.SumSquare())/2.;
//   	u2 = vel.SumSquare();

//  	if(u2<=1e-6*(2. * UU * kond)) {
//  	  Peclet=0.;
//  	  psi=0.;
//  	  tau=0.;
//  	} else {
//  	  Peclet = u2 / (2. * UU * kond);
//  	  psi = 1./tanh(Peclet)-1/Peclet;
//  	  tau = psi/(2.*UU);
//  	}
//   	Pert  = tau * (vel.Row(kdof+1)) * DSHAPEX;
//   	W = SHAPE + Pert;
 
//   	//  	if (ijob==COMP_MAT) matloc += WPG * (kond * DSHAPEX.t() * DSHAPEX
//   	//  					     + W.t() * vel * DSHAPEX);
//   	veccontr.Column(kdof+1)
//   	  +=  WPG * (kond * DSHAPEX.t()
//  		     + W.t() * (vel.Row(kdof+1)) ) * grad_fi.Column(kdof+1);
//        }
//      }
//      veccontr >> &(RETVAL(k,0,0));


      ierr = FMp(tmp,dshapext,grad_fi); FMCHK;
      ierr = FMaxpy(veccontr,WPG*kond,tmp); FMCHK;
    }
    veccontr.copy(&(RETVAL(k,0,0)));
  }
  for (ipg=0; ipg<npg; ipg++) {
    delete FM_shape[ipg];
    delete FM_dshapexi[ipg];
  }
  delete[] FM_shape;
  delete[] FM_dshapexi;

  //  exit(0);
}

#undef SHAPE    
#undef DSHAPEXI 
#undef DSHAPEX  
#undef WPG      

# Local Variables: $
# mode: c++ $
# End: $
*/
