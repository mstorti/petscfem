//__INSERT_LICENSE__
//$Id: srfgath.cpp,v 1.1 2004/01/26 20:22:34 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "./srfgath.h"

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
  n.setel(1.,1);
  n.setel(1.,2);
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
  FastMat2 Jaco(2,ndimel,ndim),staten(2,nel,ndof), 
    stateo(2,nel,ndof),u_old(1,ndof),u(1,ndof),
    n(1,ndim),xpg(1,ndim);

  Time * time_c = (Time *)time;
  double t = time_c->time();
  
  // Initialize the call back functions
  init();

  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  int ielh=-1,kloc;
  for (int k=el_start; k<=el_last; k++) {
    if (!compute_this_elem(k,this,myrank,iter_mode)) continue;
    FastMat2::reset_cache();
    ielh++;

    for (kloc=0; kloc<nel; kloc++) {
      int node = ICONE(k,kloc);
      xloc.ir(1,kloc+1).set(&NODEDATA(node-1,0));
    }
    xloc.rs();

    staten.set(&(LOCST(ielh,0,0)));
    stateo.set(&(LOCST2(ielh,0,0)));

    GPdata::edge q;
    int e=0;
    for (q = gp_data.edges_begin(); q!=gp_data.edges_end(); q++) {
      printf("edge %d, nodes %d-%d\n",e++,q.first(),q.second());
    }
    PetscFinalize();
    exit(0);
 
#if 0
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
  FastMat2::void_cache();
  FastMat2::deactivate_cache();
  return 0;
}

#undef SHAPE    
#undef DSHAPEXI 
#undef WPG      
#undef SQ
