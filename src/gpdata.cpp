//__INSERT_LICENSE__
//$Id: gpdata.cpp,v 1.38 2004/01/29 00:12:44 mstorti Exp $

#include "petscsles.h"
#include <math.h>
#include <fnmatch.h>

#include <src/fem.h>
#include <src/util2.h>
#include <src/utils.h>
#include <src/gpdata.h>

#define GPERROR \
    {PFEM_TRACE(""); \
    PetscPrintf(PETSC_COMM_WORLD,"Not implemented combination: geometry=\"%s\","\
           "nel=%d, npg=%d, ndimel=%d, \n", \
	   geom,nel,npg,ndimel); \
    PetscFinalize(); \
		       exit(0);} \

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "cartesian_2d_shape(Rowvector &,Matrix &,matrix &)"
int cartesian_2d_shape(RowVector &shape,Matrix &dshapexi,
		       double xipg,double etapg) {

  int nel=4,ndim=2;
  static int flag=0;
  static Matrix xinode;
  if (flag==0) {
    flag=1;
    xinode.ReSize(2,4);
    xinode << -1 << 1 << 1 << -1 << -1 << -1 << 1 << 1;
  }

  shape.ReSize(nel);
  dshapexi.ReSize(ndim,nel);

  for (int iloc=1; iloc<=nel; iloc++) {
    // 1D Shape Functions
    // along X
    double sxi=(1.+xipg*xinode(1,iloc))/2.;
    double dsxidxi=xinode(1,iloc)/2.;

    // along Y
    double seta=(1.+etapg*xinode(2,iloc))/2.;
    double dsetadeta=xinode(2,iloc)/2.;
	  
    shape(iloc)=sxi*seta;
    dshapexi(1,iloc)=dsxidxi*seta;
    dshapexi(2,iloc)=sxi*dsetadeta;
  }
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double xipgf(int ipg,int npg1d) {
  if (npg1d==2) {
    return (2*ipg-1)/sqrt(3.);
  } else {
    return 0.;
  }
}

void cart_prod(int npg,int nel,int nel_lay,
	       int ndim,RowVector *shape,Matrix *dshapexi,Matrix *dshapex,double *wpg,
	       GPdata &base_gp) {
  for (int ipg=0; ipg<npg; ipg++) {
    shape[ipg] = RowVector(nel);
    shape[ipg].Columns(1,nel_lay) = base_gp.shape[ipg];

    dshapexi[ipg]= Matrix(ndim,nel);
    dshapexi[ipg]= 0.;
    dshapexi[ipg].SubMatrix(1,ndim-1,1,nel_lay) = base_gp.dshapexi[ipg];
    assert(nel % nel_lay ==0);
    int nlay = nel/nel_lay;
    double coef2[] = {-0.5, 0.5};
    double coef3[] = {-1.5,2.0,-0.5};
    double coef4[] = {-(11./6.), +(18./6.), -( 9./6.), +( 2./6.)};
    double *coef;
    if (nlay == 2) coef = coef2;
    else if (nlay==3) coef = coef3;
    else if (nlay==4) coef = coef4;
    else assert(0);
#define NEWMAT_DIMS(m) printf(#m ": dims -> %d x %d\n", \
			      (m).Nrows(),(m).Ncols())
    for (int lay=0; lay<nlay; lay++)
      // NEWMAT_DIMS(dshapexi[ipg]);
      dshapexi[ipg].	
	SubMatrix(ndim,ndim,lay*nel_lay+1,(lay+1)*nel_lay) = 
	coef[lay]*base_gp.shape[ipg];
    
    wpg[ipg] = base_gp.wpg[ipg];
     dshapex[ipg]= Matrix(ndim,nel);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool operator!=(GPdata::edge p,GPdata::edge q) { 
  assert(q.gp==p.gp); 
  return p.indx!=q.indx; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GPdata::edge::edge(int j,const GPdata *gp_a) : indx(j), gp(gp_a) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GPdata::edge::edge() : indx(-1), gp(NULL) { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GPdata::edge 
GPdata::edge::operator++(int) { 
  edge q = *this;
  indx++;
  return q; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int GPdata::edge::first() { assert(gp); return gp->edges[indx*2]; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int GPdata::edge::second() { assert(gp); return gp->edges[indx*2+1]; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GPdata::edge 
GPdata::edges_begin() { return edge(0,this); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GPdata::edge 
GPdata::edges_end() { return edge(nedges_m,this); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int GPdata::nedges() { return nedges_m; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GPdata::GPdata(const char *geom,int ndimel,int nel,int npg_,int
	       mat_version_) {
  mat_version = mat_version_;
  npg= npg_;
  int ipg;
  nedges_m = 0;
  // Dimension GP vectors

  // ndimel:= dimension of the element. Cartesian elements are the
  //          cartesian product of 1D linear elements. ndimel is not
  //          necessarily equal to the number of dimensions in the
  //          space. For instance, a boundary element (for imposed
  //          flux or convection boundary condition) has a dimension
  //          lower than the spatial dimension.

  // that's me....
  // this is no true for elements that are points (bccconvs for 1D elements).

  // Shape function
  shape = new RowVector[npg];
  // Weights
  wpg = new double[npg];
  // dshapexi:= [ipg](jd,jloc) Gradient of shape function
  //    with respect to master element coordinates
  //    ipg = number of Gauss Point
  //    jd = spatial coordinate
  //    jloc = local node number
  dshapexi = new Matrix[npg];
  // dshapexi:= [ipg](jd,jloc) Gradient of shape function
  //    with respect to global coordinates
  dshapex = new Matrix[npg];

  if ( !(strcmp(geom,"prismatic")) ) {
    assert(nel==6);

    // A prisma with triangular base, useful when extruding
    // triangular meshes
    master_volume = 1;		// that is 0.5(tri)*2(1d-segment)
    assert(ndimel==3);
    assert(nel==6);
    assert(npg==1 || npg==6 || npg==8); // other cases may be considered
    int npg_seg, npg_tri;
    if (npg==1) { npg_seg=1; npg_tri=1; }
    else if (npg==6) { npg_seg=2; npg_tri=3; }
    else if (npg==8) { npg_seg=2; npg_tri=4; }
    GPdata gp_seg("cartesian1d",1,2,npg_seg);
    GPdata gp_tri("triangle",2,3,npg_tri);

    int ipg = 0;
    for (int ipg_seg=0; ipg_seg<npg_seg; ipg_seg++) {
      for (int ipg_tri=0; ipg_tri<npg_tri; ipg_tri++,ipg++) {
	// Weights are simply the poduct of the seg/tri weights
	wpg[ipg] = gp_seg.wpg[ipg_seg]*gp_tri.wpg[ipg_tri];
	// Shape fuction are the kronecker product also
	shape[ipg] = kron(gp_seg.shape[ipg_seg],gp_tri.shape[ipg_tri]);
	// x,y gradients of shape functions are the poduct of the
	// gradients of the triangle and the shape function of the segment

	dshapexi[ipg]= Matrix(ndimel,nel);
	dshapexi[ipg].Rows(1,2) = kron(gp_seg.shape[ipg_seg],gp_tri.dshapexi[ipg_tri]);
	// `z' gradients of shape functions are the poduct of the
	// shape functions of the triangle and the gradient of the
	// shape function of the segment
	dshapexi[ipg].Row(3) = kron(gp_seg.dshapexi[ipg_seg],gp_tri.shape[ipg_tri]);
	// global gradients have only to be defined (they are computed with the actual
	// coordinates. (this should be done once fo all geometries. Currently the
	// code is duplicated for each geometry :-(  )
	dshapex[ipg]= Matrix(ndimel,nel);
      }
    }
#ifdef USE_DX
    splitting.parse("1  2  3  4   5  4  6  2   2  6  3  4 tetrahedra");
#endif
    // edges
    nedges_m = 9;
    // edges.resize(2*nedges_m);
    int edges_v[] = {1,2, 2,3, 3,1, // top
		     4,5, 5,6, 6,4,
		     1,4, 2,5, 3,6};
    edges.insert(edges.end(),edges_v,edges_v+2*nedges_m);

  } else if ( !(strcmp(geom,"triangle")) ) {
    assert(nel==3);
    master_volume = 0.5;

    double xipg,etapg;
    for (ipg=0; ipg<npg; ipg++) {
      if (npg==7) {
	// Three points in the edge centers
#define A1 (0.0597158717)
#define B1 (0.4701420641)
#define A2 (0.7974269853)
#define B2 (0.1012865073)

#define W1 (0.225)
#define W2 (0.1323941527)
#define W3 (0.1259391805)

	double wpg_[7]={W1,W2,W2,W2,W3,W3,W3};
	double xpg_[7][2]={
	  1./3., 1./3.,
	  A1, B1,
	  B1, A1,
	  B1, B1,
	  A2, B2,
	  B2, A2,
	  B2, B2};
	xipg  = xpg_[ipg][0];
	etapg = xpg_[ipg][1];
	wpg[ipg] = wpg_[ipg];

      } else if (npg==4) {
	// Three points in the edge centers
	double wpg_[4]={-27./96.,25./96.,25./96.,25./96.};
	double xpg_[4][2]={
	  1./3., 1./3.,
	  0.2, 0.2,
	  0.6, 0.2,
	  0.2, 0.6};
	xipg  = xpg_[ipg][0];
	etapg = xpg_[ipg][1];
	wpg[ipg] = wpg_[ipg];

      } else if (npg==3) {

	// Three points in the edge centers
	double a=0.5;
	if (ipg==0) {
	  xipg=1-2*a; etapg=a;
	} else if (ipg==1) {
	  xipg=a; etapg=1-2*a;
	} else {
	  xipg=a; etapg=a;
	}
	wpg[ipg] =  1/6.;

      } else if (npg==1) {
	// One point in the center of the element
	xipg=1./3.;
	etapg=1./3.;
	wpg[ipg] =  0.5;

      } else {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Not valid value of npg. for triangles %d\n",
		    npg);
	PetscFinalize();
	exit(0);
      }

      shape[ipg] = RowVector(nel);
      shape[ipg](1)=xipg;
      shape[ipg](2)=etapg;
      shape[ipg](3)=1-xipg-etapg;
    
      dshapex[ipg]= Matrix(ndimel,nel);
      dshapexi[ipg]= Matrix(ndimel,nel);
    
      dshapexi[ipg](1,1)=1.;
      dshapexi[ipg](2,1)=0.;
      dshapexi[ipg](1,2)=0.;
      dshapexi[ipg](2,2)=1.;
      dshapexi[ipg](1,3)=-1.;
      dshapexi[ipg](2,3)=-1.;
    }

#ifdef USE_DX
    splitting.parse("1  2  3  triangles");
#endif
  } else if ( !(strcmp(geom,"tetra")) ) {
    assert(nel==4);
    master_volume = 1./6.;
    double xipg,etapg,zetapg;
    for (ipg=0; ipg<npg; ipg++) {
      if (npg==4) {
	// Four points near  the side centers
	double a=(30.+sqrt(660.))/120.;
	xipg = a; etapg=a; zetapg=a; 
	if (ipg==0) {
	  xipg=1-3*a; 
	} else if (ipg==1) {
	  etapg=1-2*a;
	} else {
	  zetapg=1-2*a;
	}
	wpg[ipg] =  master_volume/4.0;

      } else if (npg==1) {
	// One point in the center of the element
	xipg=1./3.;
	etapg=1./3.;
	zetapg=1./3.;
	wpg[ipg] =  master_volume;

      } else {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Not valid value of npg. for tetras %d\n",
		    npg);
	PetscFinalize();
	exit(0);
      }

      shape[ipg] = RowVector(nel);
      shape[ipg](1)=1-xipg-etapg-zetapg;
      shape[ipg](2)=xipg;
      shape[ipg](3)=etapg;
      shape[ipg](4)=zetapg;
    
      dshapex[ipg]= Matrix(ndimel,nel);
      dshapexi[ipg]= Matrix(ndimel,nel);
    
      dshapexi[ipg]=0;
      dshapexi[ipg](1,2)=1.;
      dshapexi[ipg](2,3)=1.;
      dshapexi[ipg](3,4)=1.;
      dshapexi[ipg](1,1)=-1.;
      dshapexi[ipg](2,1)=-1.;
      dshapexi[ipg](3,1)=-1.;
    }
    // edges
    nedges_m = 6;
    // edges.resize(2*nedges_m);
    int edges_v[] = {1,2, 2,3, 3,1, 4,1, 4,2, 4,3};
    edges.insert(edges.end(),edges_v,edges_v+2*nedges_m);

#ifdef USE_DX
    splitting.parse("1  2  3  4  tetrahedra");
#endif

  } else if (!fnmatch("cartesian*d_face",geom,0)) {

    // This is for integration over the face of an
    // element. (For instance in BCCONV elements.)

    // npg1d:= number of points per dimensional direction
    // npg:= npg1d^ndimel total number of Gauss points
    int ndimel;
    sscanf(geom,"cartesian%dd",&ndimel);
    master_volume = pow(2.0,ndimel);

//      if (ndimel==0 && npg==1) {
    int npg1d=int(pow(double(npg),1./double(ndimel)));
    // AGREGAR lin1d y brick integrations!!
    //    } else if (!(strcmp(geom,"quad")) || !(strcmp(geom,"lin1d"))
    //  	     || !(strcmp(geom,"brick")) ) {

    if (npg1d!=2) GPERROR;
    // Position of the nodes in the master element
    
    if (ndimel==2) {

      ipg=-1;
      for (int ixipg=0; ixipg<=1; ixipg++) {
	double xipg=(2*ixipg-1)/sqrt(3.);
	double etapg = -1;
	ipg++;
	wpg[ipg] = 1;

	cartesian_2d_shape(shape[ipg],dshapexi[ipg],xipg,etapg);
      }
    } else { GPERROR; } 
#ifdef USE_DX
    splitting.parse("1  2  4  3 quads");
#endif
  } else if (!fnmatch("cartesian*d",geom,0)) {

    // npg1d:= number of points per dimensional direction
    // npg:= npg1d^ndimel total number of Gauss points
    int ndimel;
    sscanf(geom,"cartesian%dd",&ndimel);
    assert(nel==int_pow(2,ndimel));
    master_volume = pow(2.0,ndimel);
    int npg1d=int(pow(double(npg),1./double(ndimel)));
    // AGREGAR lin1d y brick integrations!!
    //    } else if (!(strcmp(geom,"quad")) || !(strcmp(geom,"lin1d"))
    //  	     || !(strcmp(geom,"brick")) ) {

    if (npg1d!=2 && npg1d!=1 && ndimel>0) GPERROR;
    
    if (ndimel==0 && npg==1) {
      // for bccconv 0D elements npg=1, ndimel=0, nel=1, geom=cartesian0d.
      // the value for shape function is a ndim x nel x npg vector equal to 1.0 at this point 
      // and gradient of shape func is undefined (then dim=0)
      
      shape[0] = RowVector(nel);
      shape[0](1)=1.;
      wpg[0]=1.;
      dshapex[0] = Matrix(ndimel,nel);
      dshapexi[0] = Matrix(ndimel,nel);

    } else if (ndimel==1) {
      // Position of the nodes in the master element
      Matrix xinode(1,2);
      xinode << -1 << 1 ;
      
      ipg=-1;
      for (int ixipg=0; ixipg<npg1d; ixipg++) {
	double xipg = xipgf(ixipg,npg1d);
	ipg++;
	wpg[ipg] = master_volume/double(npg);
	
	shape[ipg] = RowVector(nel);
	dshapex[ipg]= Matrix(ndimel,nel);
	dshapexi[ipg]= Matrix(ndimel,nel);
	
	for (int iloc=1; iloc<=nel; iloc++) {
	  // 1D Shape Functions
	  // along X
	  double sxi=(1.+xipg*xinode(1,iloc))/2.;
	  double dsxidxi=xinode(1,iloc)/2.;

	  shape[ipg](iloc)=sxi;
	  dshapexi[ipg](1,iloc)=dsxidxi;
	}
      }
    } else if (ndimel==2) {
      Matrix xinode(2,4);
      xinode << -1 << 1 << 1 << -1 << -1 << -1 << 1 << 1;
      
      ipg=-1;
      for (int ixipg=0; ixipg<npg1d; ixipg++) {
	double xipg = xipgf(ixipg,npg1d);
	for (int ietapg=0; ietapg<npg1d; ietapg++) {
	  ipg++;
	  double etapg=xipgf(ietapg,npg1d);
	  wpg[ipg] = master_volume/double(npg);

	  shape[ipg] = RowVector(nel);
	  dshapex[ipg]= Matrix(ndimel,nel);
	  dshapexi[ipg]= Matrix(ndimel,nel);
	
	  for (int iloc=1; iloc<=nel; iloc++) {
	    // 1D Shape Functions
	    // along X
	    double sxi=(1.+xipg*xinode(1,iloc))/2.;
	    double dsxidxi=xinode(1,iloc)/2.;

	    // along Y
	    double seta=(1.+etapg*xinode(2,iloc))/2.;
	    double dsetadeta=xinode(2,iloc)/2.;

	    shape[ipg](iloc)=sxi*seta;
	    dshapexi[ipg](1,iloc)=dsxidxi*seta;
	    dshapexi[ipg](2,iloc)=sxi*dsetadeta;
	  }
	}
      }
#ifdef USE_DX
    splitting.parse("1  2  4  3 quads");
#endif
    } else if (ndimel==3) {
      Matrix xinode(3,8);
      xinode << -1 << 1 << 1 << -1 << -1 << 1 << 1 << -1 
	     << -1 << -1 << 1 << 1 << -1 << -1 << 1 << 1 
	     << -1 << -1 << -1 << -1 << 1 << 1 << 1 << 1;
      
      ipg=-1;
      for (int ixipg=0; ixipg<npg1d; ixipg++) {
	double xipg = xipgf(ixipg,npg1d);
	for (int ietapg=0; ietapg<npg1d; ietapg++) {
	  double etapg=xipgf(ietapg,npg1d);
	  for (int izetapg=0; izetapg<npg1d; izetapg++) {
	    double zetapg=xipgf(izetapg,npg1d);

	    ipg++;
	    wpg[ipg] = master_volume/double(npg);

	    shape[ipg] = RowVector(nel);
	    dshapex[ipg]= Matrix(ndimel,nel);
	    dshapexi[ipg]= Matrix(ndimel,nel);
	    
	    for (int iloc=1; iloc<=nel; iloc++) {
	      // 1D Shape Functions
	      // along X
	      double sxi=(1.+xipg*xinode(1,iloc))/2.;
	      double dsxidxi=xinode(1,iloc)/2.;
	      
	      // along Y
	      double seta=(1.+etapg*xinode(2,iloc))/2.;
	      double dsetadeta=xinode(2,iloc)/2.;
	      
	      // along Z
	      double szeta=(1.+zetapg*xinode(3,iloc))/2.;
	      double dszetadzeta=xinode(3,iloc)/2.;
	      
	      shape[ipg](iloc)=sxi*seta*szeta;
	      dshapexi[ipg](1,iloc)=dsxidxi*seta*szeta;
	      dshapexi[ipg](2,iloc)=sxi*dsetadeta*szeta;
	      dshapexi[ipg](3,iloc)=sxi*seta*dszetadzeta;
	    }
	  }
	}
      }
#ifdef USE_DX
    splitting.parse("1 2 4 3  5 6 8 7 cubes");
#endif
    // edges
    nedges_m = 12;
    int edges_v[] = {1,2, 2,3, 3,4, 4,1, // Bottom
		     5,6, 6,7, 7,8, 8,5, //Top
		     1,5, 2,6, 3,7, 4,8 }; // Lateral
    edges.insert(edges.end(),edges_v,edges_v+2*nedges_m);

    } else GPERROR;
      
  } else if (!strcmp(geom,"line2quad")) {
    
    assert(ndimel==2);
    assert(nel % 2==0);
    int nlay = nel/2;
    assert(nlay>=2);
    assert(nlay<=4);
    assert(npg==2); // other cases may be considered
    GPdata line("cartesian1d",1,2,npg,GP_NEWMAT);
    cart_prod(npg,nel,2,ndimel,shape,dshapexi,dshapex,wpg,line);

  } else if (!strcmp(geom,"quad2hexa")) {
    
    assert(ndimel==3);
    assert(nel % 4 == 0 && 2<=(nel/4) && (nel/4)<=4 );
    assert(npg==4); // other cases may be considered
    GPdata quad("cartesian2d",ndimel-1,4,npg,GP_NEWMAT);
    cart_prod(npg,nel,4,ndimel,shape,dshapexi,dshapex,wpg,quad);

  } else if (!strcmp(geom,"tri2prism")) {
    
    assert(ndimel==3);
    assert(nel % 3 == 0 && 2<=(nel/3) && (nel/3)<=4 );
    assert(npg==3); // other cases may be considered
    GPdata tri("triangle",ndimel-1,3,npg,GP_NEWMAT);
    cart_prod(npg,nel,3,ndimel,shape,dshapexi,dshapex,wpg,tri);

  } else GPERROR;
  
  if (mat_version==GP_FASTMAT) {

    // make fastmat copy version
    FM_shape = new (FastMat *)[npg];
    FM_dshapexi = new (FastMat *)[npg];
    for (int ipg=0; ipg<npg; ipg++) {
      FM_shape[ipg] = new FastMat;
      NM2FM(*FM_shape[ipg],shape[ipg]);
      FM_dshapexi[ipg] = new FastMat;
      NM2FM(*FM_dshapexi[ipg],dshapexi[ipg]);
    }

  } else if (mat_version==GP_FASTMAT2) {

    FM2_shape = new (FastMat2 *)[npg];
    FM2_dshapexi = new (FastMat2 *)[npg];
    for (int ipg=0; ipg<npg; ipg++) {
      FM2_shape[ipg] = new FastMat2;
      FM2_shape[ipg]->set(shape[ipg]);
      int nel = shape[ipg].Ncols();
      int nrows = shape[ipg].Nrows();
      assert(nrows==1);
      FM2_shape[ipg]->reshape(1,nel);
      FM2_dshapexi[ipg] = new FastMat2;
      FM2_dshapexi[ipg]->set(dshapexi[ipg]);
    }

  }
  wpg_sum = 0.;
       for (int ipg=0; ipg<npg; ipg++) {
	 wpg_sum += wpg[ipg];
       }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GPdata::~GPdata() {
  delete[] wpg;
  delete[] shape;
  delete[] dshapexi;
  delete[] dshapex;
  if (mat_version == GP_FASTMAT) {
    for (int ipg=0; ipg<npg; ipg++) {
      delete FM_shape[ipg];
      delete FM_dshapexi[ipg];
      // FM_shape[ipg]->~FastMat();
      // FM_dshapexi[ipg]->~FastMat();
    }
    delete[] FM_shape;
    delete[] FM_dshapexi;
  } else if (mat_version == GP_FASTMAT2) {
    for (int ipg=0; ipg<npg; ipg++) {
      delete FM2_shape[ipg];
      delete FM2_dshapexi[ipg];
      // FM_shape[ipg]->~FastMat();
      // FM_dshapexi[ipg]->~FastMat();
    }
    delete[] FM2_shape;
    delete[] FM2_dshapexi;
  }
}

#undef GPERROR
