//__INSERT_LICENSE__
//$Id: srfgath.cpp,v 1.3 2004/01/27 19:16:42 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "./srfgath.h"

enum sf_error { not_correct_number_of_intersections };
SurfGatherer::SurfFunction::~SurfFunction() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class plane : public SurfGatherer::SurfFunction {
private:
  int ndim;
  FastMat2 x0,n,dx,tmp;
public:
  void init(const TextHashTable *thash);
  double f(const FastMat2 &x);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void plane::init(const TextHashTable *thash) {
  ndim = 3;
  x0.resize(1,ndim);
  n.resize(1,ndim);
  dx.resize(1,ndim);
  x0.set(0.);
  n.set(0.).setel(1.,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double plane::f(const FastMat2 &x) {
  dx.set(x).rest(x0);
  tmp.prod(dx,n,-1,-1);
  return double(tmp);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
SurfGatherer::SurfFunction* 
SurfGatherer::SurfFunction::factory(const TextHashTable *thash) {
  int ierr;
  SurfGatherer::SurfFunction *sf;
  //o Defines the geomtry of the element
  TGETOPTDEF_S(thash,string,surf_fun_type,<none>);
  assert(surf_fun_type!="<none>");
  if (surf_fun_type =="plane") {
    sf = new plane;
  } else 
    PETSCFEM_ERROR("SurfGatherer::SurfFunction::factory: "
		   "unknown surf_fun_type \"%s\"\n",surf_fun_type.c_str());
  return sf;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void SurfGatherer::initialize() {
  int ierr;
  sf = SurfGatherer::SurfFunction::factory(thash);
  sf->init(thash);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
SurfGatherer::~SurfGatherer() { if (sf) delete sf; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int SurfGatherer::ask"
int SurfGatherer::ask(const char *jobinfo,int &skip_elemset) {
  skip_elemset = 1;
  DONT_SKIP_JOBINFO(gather);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "SurfGatherer::assemble"
int SurfGatherer::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			   Dofmap *dofmap,const char *jobinfo,int myrank,
			   int el_start,int el_last,int iter_mode,
			   const TimeData *time) {
  
  int ierr;

  GET_JOBINFO_FLAG(gather);
  assert(gather);

  //o Dimension of embedding space
  TGETOPTDEF(thash,int,ndim,0);
  assert(ndim==3);
  //o Position in gather vector
  TGETOPTDEF(thash,int,gather_pos,0);
  //o How many gather values will be contributed by this elemset per surface. 
  TGETOPTDEF(thash,int,vals_per_plane,0);
  assert(vals_per_plane>0);
  //o Number of Gauss points.
  TGETOPTNDEF(thash,int,npg,none);
  //o Defines the geomtry of the element
  TGETOPTDEF_S(thash,string,geometry,<none>);
  assert(geometry!="<none>");

  //o Values #c_j# where #f(x,y,z)=c_j# surfaces are defined. 
  // _T: double array
  // _N: f_vals
  // _D: 0.0
  // _DOC: If an array #c_0# , #c_1# ,.... #c_{n-1}# is entered, then
  // the quantities are integrated on surfaces for #f(x,y,z)=c_0# ,
  // #f(x)=c_1# , ..., #f(x)=c_{n-1}# . 
  // _END
  vector<double> f_vals;
  const char *line;
  thash->get_entry("f_vals",line);
  if (line) {
    read_double_array(f_vals,line);
    assert(f_vals.size()>0);
  } else f_vals.resize(1,0);
  int nvals = f_vals.size();
  int gather_length = vals_per_plane*nvals;
  
  GPdata gp_data(geometry.c_str(),ndim,nel,npg,GP_FASTMAT2);

#define NODEDATA(j,k) VEC2(nodedata->nodedata,j,k,nu)
#define ICONE(j,k) (icone[nel*(j)+(k)]) 

#define LOCST(iele,j,k) VEC3(locst,iele,j,nel,k,ndof)
#define LOCST2(iele,j,k) VEC3(locst2,iele,j,nel,k,ndof)

  // get number of fields per node (constant+variables)
  int nu=nodedata->nu;

  int ja = 0;
  double *locst = arg_data_v[ja++].locst;
  double *locst2 = arg_data_v[ja++].locst;
  int options = arg_data_v[ja].options;
  vector<double> *values = arg_data_v[ja++].vector_assoc;
  int nvalues = values->size();
  // check that we don't put values beyond the end of global vector
  // `values'
  assert(gather_pos + gather_length <= nvalues); 

  FastMat2 xloc(2,nel,ndim);

  int ndimel = ndim-1;
  int nedges = gp_data.nedges();
  assert(nedges>0);
  FastMat2 Jaco(2,ndimel,ndim),staten(2,nel,ndof), 
    stateo(2,nel,ndof),u_old(1,ndof),u(1,ndof),
    n(1,ndim),xpg(1,ndim),x(1,ndim),xi(2,ndim,nedges),
    ui(2,ndof,nedges), xc(1,ndim), xcp(1,ndim), uc(1,ndof),
    x1(1,ndim), x2(1,ndim), dx(1,ndim), tmp;
  xi.set(0.);

  Time * time_c = (Time *)time;
  double t = time_c->time();
  
  // Initialize the call back functions
  init();

  FastMatCacheList cache_list;
  // FastMat2::activate_cache(&cache_list);
  vector<double> f(nel), alpha(nedges);
  vector<int> indx(nedges);

  int ielh=-1,kloc;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;

    double max_f, min_f;
    for (kloc=0; kloc<nel; kloc++) {
      int node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
      x.set(xloc);
      f[kloc] = sf->f(x);
      if (kloc==0) {
	max_f = f[0];
	min_f = f[0];
      }
      if (f[kloc]>max_f) max_f = f[kloc];
      if (f[kloc]<min_f) min_f = f[kloc];
    }
    xloc.rs();

    for (int jval=0; jval<nvals; jval++) {
      // Compute surface `jval', given by `f=val'
      double val = f_vals[jval];
      // Check if surface intersects element
      if (max_f<val || min_f>val) continue;
      // Nbr of intersections (and number of sides
      // of the polygon)
      int nint = 0;
      for (GPdata::edge q = gp_data.edges_begin(); 
	   q!=gp_data.edges_end(); q++) {
	// compute intersection with edge
	int n1=q.first();
	int n2=q.second();
	double f1 = f[n1-1];
	double f2 = f[n2-1];
	// check if surface intersects edge 
	if ((f1-val)*(f2-val)>=0.0) continue;
	int jint = nint++;
	// intersection is at (1-beta)*x1+beta*x1
	double beta = (val-f1)/(f2-f1); 
	// xi(:,j) := coordinates of intersection j
	xi.ir(2,jint+1);
	xloc.ir(1,n1);
	xi.set(xloc).scale(1-beta);
	xloc.ir(1,n2);
	xi.axpy(xloc,beta);
	// ui(:,j) := state at intersection j
	ui.ir(2,jint+1);
	staten.ir(1,n1);
	ui.set(staten).scale(1-beta);
	staten.ir(1,n2);
	ui.axpy(staten,beta);
      }
      xi.rs();
      ui.rs();
#if 0
      printf("elem %d, intersections: %d\n",k,nint);
      if(nint>0) {
	xi.is(2,1,nint).print("intersections: ");
	xi.rs();
      }
#endif
      if (nint<3) set_error(not_correct_number_of_intersections);
      //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
      // Intersection area may be a triangle or quad for a tetra. Or
      // higher polygons in other cases. We assume that it is a
      // convex polygon and integrate over them by dividing the
      // polygon in triangles from the center of the polygon to
      // each side. 

      // Compute center of polygon xc
      xi.is(2,1,nint).print("intersections: ");
      xc.sum(xi,1,-1).scale(1./double(nint));
      // Computes normal as gradient of f in the center
      double epsil = 1e-3,fp,fm;
      for (int j=1; j<=ndim; j++) {
	xcp.set(xc).addel(epsil,j);
	fp = sf->f(xcp);
	xcp.addel(-2*epsil,j);
	fm = sf->f(xcp);
	n.setel((fp-fm)/(2*epsil),j);
      }
      double an = n.norm_p_all(2.0);
      n.scale(1./an);

      // uc:= value of u interpolated at center of polygon
      ui.is(2,1,nint);
      uc.sum(ui,1,-1).scale(1./double(nint));
      // x1,x2 perpendicular to normal n
      // are two vectors that define the axis on the
      // plane normal to the polygon
      xi.rs().ir(2,1);
      x1.cross(n,xi);
      double l1 = x1.norm_p_all(2.0);
      x1.scale(1./l1);
      x2.cross(n,x1);
      // For each vertex of the polygon, we compute an angle
      // of rotation around it, in order to project. 
      for (int j=0; j<nint; j++) {
	xi.ir(2,j+1);
	dx.set(xi).rest(xc);
	tmp.prod(x1,dx,-1,-1);
	double X1 = tmp.get();
	tmp.prod(x2,dx,-1,-1);
	double X2 = tmp.get();
	alpha[j] = atan2(X2,X1);
      }
      xi.rs();

      for (int j=0; j<nint; j++) indx[j] = j;
#if 1
      printf("sorting angles\n");
      for (int j=0; j<nint; j++) printf("(%d,%g) ",j,alpha[j]);
      printf("\n");
#endif
      // Sort angles and determine permutation
      for (int j=0; j<nint-1; j++) {
	// Look for the minimum in range j,nint and exchange with pos. j
	double amin = alpha[j];
	int kmin = j;
	for (int k=j+1; k<nint; k++) {
	  if (alpha[k]<amin) {
	    amin = alpha[k];
	    kmin = k;
	  }
	}
	alpha[kmin] = alpha[j];
	alpha[j] = amin;
	indx[kmin] = j;
	indx[j] = kmin;
      }
#if 1
      for (int j=0; j<nint; j++) printf("(%d,%g) ",indx[j],alpha[j]);
      printf("\n");
#endif
    }
 
#if 0
    staten.set(&(LOCST(ielh,0,0)));
    stateo.set(&(LOCST2(ielh,0,0)));

    GPdata::edge q;
    int e=0;
    
    for (q = gp_data.edges_begin(); q!=gp_data.edges_end(); q++) {
      int n1=q.first();
      int n2=q.second();
      double f1 
    }
 
    // Let user do some things when starting with an element
    element_hook(k);

    for (int ipg=0; ipg<npg; ipg++) {
      // Gauss point coordinates
      xpg.prod(SHAPE,xloc,-1,-1,1);
      // Jacobian master coordinates -> real coordinates
      Jaco.prod(DSHAPEXI,xloc,1,-1,-1,2);

      double detJaco;
      if (ndimel==ndim) {
	detJaco = Jaco.det();
      } else if (ndimel==ndim-1) {
	detJaco = Jaco.detsur(&n);
	n.scale(1./detJaco);
	n.scale(-1.);		// fixme:= This is to compensate a bug in mydetsur
      }
      if (detJaco <= 0.) {
	printf("Jacobian of element %d is negative or null\n"
	       " Jacobian: %f\n",k,detJaco);
	PetscFinalize();
	exit(0);
      }
      double wpgdet = detJaco*WPG;

      // Values of variables at Gauss point
      u.prod(SHAPE,staten,-1,-1,1);
      u_old.prod(SHAPE,stateo,-1,-1,1);

      set_pg_values(pg_values,u,u_old,xpg,n,wpgdet,t);
      if (options & VECTOR_ADD) {
	for (int j=0; j<gather_length; j++) {
	  (*values)[gather_pos+j] += pg_values[j];
	}
      } else assert(0);
    }
#endif
  }  
  PetscFinalize();
  exit(0);
 
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

void SurfGatherer::handle_error(int error) {
  string s = "Unknown error";
  if (error==not_correct_number_of_intersections) 
    s = "Not correct number of intersection";
  PETSCFEM_ERROR("Elemset name \"%s\", ptr %p, set error %d\n",
		 name(),this,error);  
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
