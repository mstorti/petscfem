//__INSERT_LICENSE__
//$Id: interpola.cpp,v 1.3 2003/06/09 02:37:18 mstorti Exp $

#include <src/fem.h>
#include <src/fastmat2.h>
#include <src/interpola.h>

// extern Mesh *GLOBAL_MESH;
// TODO
// * set `profile' to 0

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "interpolation::read"
void interpolation::
read(FileStack *fstack) {

  int ierr;
  //o The space dimension
  NSGETOPTDEF_ND(int,ndim,0);
  assert(ndim>0);

  // nodedata = GLOBAL_MESH->nodedata;
  char *line;
  vector<double> row;
  npoints = 0;
  while(!fstack->get_line(line)) {
    if (!strcmp(line,"__END_DATA__")) break;
    read_double_array(row,line);
    assert(row.size()==ndim);
    for (int j=0; j<ndim; j++) points.push(row[j]);
    npoints++;
  }
  // int &nnod = nodedata->nnod;
  points.defrag().reshape(2,npoints,ndim);
  point2elem.resize(npoints);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "int interpolation::ask(char *,int &)"
int interpolation::ask(const char *jobinfo,int &skip_elemset) {
   skip_elemset = !strcmp(jobinfo,jobinfo_stage());
   return 0;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "interpolation::assemble"
void interpolation::
new_assemble(arg_data_list &arg_data_v,const Nodedata *nodedata,
	     const Dofmap *dofmap,const char *jobinfo,
	     const ElementList &elemlist,
	     const TimeData *time_data) {

  int ierr,nelprops,nel,ndof;
  elem_params(nel,ndof,nelprops);

  //o The space dimension
  NSGETOPTDEF_ND(int,ndim,0);
  //o The geometry of the elemset
  NSGETOPTDEF(string,geometry,"cartesian2d");
  assert(geometry=="tetra");

  double tol=1e-10;

  int nH = nodedata->nu-ndim;
  FastMat2 xloc(2,nel,ndim), A(2,ndim,ndim), iA(2,ndim,ndim),
    x0(1,ndim), xp(1,ndim), L(1,ndim+1);
  FMatrix Hloc(nel,nH);

  FastMatCachePosition cp;
  FastMatCacheList cache_list;
  FastMat2::activate_cache(&cache_list);

  if (!initialized) {
    interp_coef.a_resize(2,npoints,nel).defrag();
    // set the flag for all points to 0
    for (int j=0; j<npoints; j++) point2elem.e(j)=-1;
    // for each element search for all points that
    // are inside it
    for (ElementIterator element = elemlist.begin(); 
	 element!=elemlist.end(); element++) {
      
      FastMat2::reset_cache();
      element.node_data(nodedata,xloc.storage_begin(),
			Hloc.storage_begin());

      xloc.ir(1,1);
      x0.set(xloc);
      for (int j=1; j<=ndim; j++) {
	xloc.ir(1,j+1);
	A.ir(2,j).set(xloc).rest(x0);
      }
      xloc.rs();
      A.rs();
      iA.inv(A);

      FastMat2::branch();
      FastMat2::choose(0);
      FastMat2::get_cache_position(cp);
      for (int p=0; p<npoints; p++) {
	if (point2elem.e(p)>=0) continue;
	FastMat2::jump_to(cp);
	L.is(1,2,ndim+1).prod(iA,xp,1,-1,-1);
	double l = L.sum_all();
	L.rs().setel(1.0-l,1);
	if (L.max_all()<=1.+tol && L.min_all()>-tol) {
	  L.export_vals(&interp_coef.e(p,0));
	  int pos_in_chunk;
	  element.position(point2elem.e(p), pos_in_chunk);
	}
      }
      FastMat2::resync_was_cached();
      FastMat2::leave();
      
      FastMat2::void_cache();
      FastMat2::deactivate_cache();
    }
    initialized = 1;
  }

}
