//__INSERT_LICENSE__
// $Id: nullvort.cpp,v 1.6 2003/03/10 20:09:38 mstorti Exp $

#include <applications/ns/nullvort.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/surf2vol.h>
#include <src/util2.h>
#include <src/linkgraph.h>
#include <src/cloud2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
null_vort_bo::null_vort_bo() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
null_vort_bo::~null_vort_bo() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void null_vort_bo::read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap) {
  /// The options table
  TextHashTable * thash = new TextHashTable;

  // Read options from data file
  thash->read(fstack);
  int ierr;

  // get options. Other options are processed in Surf2vol::factory
  //o The number of nodes per skin panel. 
  TGETOPTNDEF(thash,int,nel_surf,0);
  //o The elemset that is on the fluid side. 
  TGETOPTDEF_S(thash,string,volume_elemset,<none>);
  assert(volume_elemset!="<none>");
  //o The field to be used on the fictitious node as
  //  Lagrange multiplier
  TGETOPTNDEF(thash,int,fic_dof,<none>);
  // o The number of nodes per skin panel. 
  //TGETOPTDEF_S(thash,string,name,<none>);

  // Call the `factory' (constructor) for the Surf2Vol object. 
  // layers:= number of *element* layers, so that the number
  // of node layers is layers+1
  int identify_volume_elements, layers,
    use_exterior_normal,  ndimel;
  Surf2Vol *sv_gp_data=NULL;
  Elemset *vol_elem;
  Surf2Vol::factory(thash, volume_elemset, 
		    nel_surf, sv_gp_data, 
		    vol_elem, identify_volume_elements, layers,
		    use_exterior_normal,  ndimel);
  int nel = nel_surf*(layers+1);

  dvector<int> icone;
  icone.reshape(2,0,nel);
  vector<string> tokens;
  vector<int> row;
  map<int,int> coupl2fic;
  int nelem=0;
  char *line;
  while(!fstack->get_line(line)) {
    if (!strcmp(line,"__END_DATA__")) break;
    row.clear();
    read_int_array(row,line);
    assert(row.size()==nel_surf);
    for (int j=0; j<nel_surf; j++) icone.push(row[j]);
    for (int j=0; j<nel_surf*layers; j++) icone.push(1);
    nelem++;
  }
  while(!fstack->get_line(line)) {
    if (!strcmp(line,"__END_DATA__")) break;
    row.clear();
    read_int_array(row,line);
    assert(row.size()==2);
    coupl2fic[row[0]-1] = row[1]-1;
  }
  map<int,int>::iterator r,re;
  re = coupl2fic.end();

  icone.defrag();
  assert(nelem*nel == icone.size());
  assert(volume_elemset!="<none>");

  identify_volume_elements_fun(dofmap->nnod, nel_surf, layers,
			       nelem, icone.buff(), nel, vol_elem->nel,
			       vol_elem, sv_gp_data);
  delete sv_gp_data;
  sv_gp_data = NULL;

  // We use here the graph as a basic map<int, set<int> > object
  // to store the set of nodes neighbor to a node on the surface
  LinkGraph graph;
  graph.init(dofmap->nnod);
  // Nodes on the coupling surface (0-based)
  set<int> coupling_nodes;
  for (int j=0; j<nelem; j++) {
    for (int k=0; k<nel_surf; k++) {
      int nodek=icone.e(j,k)-1;
      coupling_nodes.insert(nodek);
      for (int l=0; l<nel_surf; l++) {
	// Store in graph (0 based)
	int nodel=icone.e(j,l)-1;
	// printf("adding to graph: %d, %d\n",nodek,nodel);
	graph.add(nodek,nodel);
	graph.add(nodel,nodek);
      }
    }
  }

#if 0
  GSet ngb;
  for (int j=0; j<dofmap->nnod; j++) {
    ngb.clear();
    graph.set_ngbrs(j,ngb);
    if (ngb.size()>0) {
      printf("%d -> (",j+1);
      GSet::iterator q, qe=ngb.end();
      for (q=ngb.begin(); q!=qe; q++) printf("%d ",*q+1);
      printf(")\n");
    }
  }
#endif

  // Numbe of nodes in th coupling surface
  int n_coupling_nodes = coupling_nodes.size();

  // coupling_nodes(surface_node,:) = nodes in the row of
  // this surface node (global, 0-based)
  dvector<int> coupling_nodes_table;
  coupling_nodes_table.a_resize(2,n_coupling_nodes,layers+1);

  // Maps a surface node to an index in coupling_nodes_table
  map<int,int> coupling_nodes_map;

  set<int>::iterator q, qe=coupling_nodes.end();
  int cn_indx=0;

  // Fill with -1 in order to know which are already filled
  for (int j=0; j<n_coupling_nodes*(layers+1); j++)
    coupling_nodes_table.ref(j) = -1;

  // Fill first layer and index map
  for (q=coupling_nodes.begin(); q!=qe; q++) {
    coupling_nodes_table.e(cn_indx,0) = *q;
    coupling_nodes_map[*q] = cn_indx;
    cn_indx++;
  }

  // Fills the remaining layers ( 1 <= layer <= layers) of the table
  icone.reshape(3,nelem,layers+1,nel_surf);
  for (int e=0; e<nelem; e++) { //loop over surface elements
    for (int sn=0; sn<nel_surf; sn++) { // loop over surface nodes
      int sf_node = icone.e(e,0,sn)-1; // a surface node
      // Check that the node is in the map
      assert(coupling_nodes_map.find(sf_node)
	     !=coupling_nodes_map.end());
      // Index in the `coupling_nodes_table'
      cn_indx = coupling_nodes_map[sf_node];
      // For each node in the row: if already loaded
      // then check that coincide, else load. 
      for (int l=1; l<=layers; l++) {
	int node = icone.e(e,l,sn)-1;
	if (coupling_nodes_table.e(cn_indx,l)==-1) 
	  coupling_nodes_table.e(cn_indx,l) = node; // not loaded
	else assert(coupling_nodes_table.e(cn_indx,l)==node); // already loaded
      }
    }
  }

#if 0 // Prints table of node rows. 
  for (cn_indx=0; cn_indx<n_coupling_nodes; cn_indx++) {
    printf("row %d: ",cn_indx);
    for (int l=0; l<=layers; l++) 
      printf("%d ",coupling_nodes_table.e(cn_indx,l));
    printf("\n");
  }
#endif

  // Number of nodes in the stencil (tangential to the surface)
  int n_stencil = 3;
  int nx = n_stencil * (layers+1);
  // This will be the connectivity table of the elemset
  dvector<int> icone_stencil;
  // The number of nodes for the Lagrange multiplier elemset
  // is the number of nodes in the stencil +1 (the fictitious node).
  // 1-based (because we have to send it to the elemset)
  icone_stencil.reshape(2,0,nx+1);
  // Neighbors (on surface) of node
  GSet ngb;
  // Constraint to be generated
  Constraint constraint;
  // List of nodes in stencil (1-based)
  dvector<int> stencil;
  // Nuumber of lagrange multipliers to be imposed
  int nelem_nlr=0;
  // Loop over nodes on the coupling surface
  for (cn_indx=0; cn_indx<n_coupling_nodes; cn_indx++) {
    // printf("computing constraint ...\n");
    // Node number (0-based)
    int node = coupling_nodes_table.e(cn_indx,0);
    // get neighbors (on the surface)
    ngb.clear();
    graph.set_ngbrs(node,ngb);
    // Check number of neighbors. We have a grid of ngb.size() nodes
    // on the surface. For each node on the surface there is a row of
    // (layers+1) nodes (including the surface node) in the direction
    // nomal to the surface. This makes a total amount of `n_cloud =
    // ngb.size()*(layers+1)' nodes.  If the number of neighbors is
    // less than 3 (in each tangent direction) we can't make a second
    // order precision approximation. So we skip the node. Give a
    // warning.
    if (ngb.size()!=n_stencil) {
      printf("No 3 ngbrs node: %d\n",node);
      continue;
    }
    // We put first the nodes on the surface
    // (the first layerxs)
    GSet::iterator q, qe=ngb.end();
    icone_stencil.push(node+1);
    for (q=ngb.begin(); q!=qe; q++) {
      // surface node
      int sf_node = *q;
      if (sf_node==node) continue;
      assert(coupling_nodes_map.find(sf_node)
	     !=coupling_nodes_map.end());
      icone_stencil.push(sf_node+1);
    }
    for (int j=0; j<n_stencil; j++) {
      int sf_node = icone_stencil.e(nelem_nlr,j)-1;
      assert(coupling_nodes_map.find(sf_node)
	     !=coupling_nodes_map.end());
      int cn_indx2 = coupling_nodes_map[sf_node];
      for (int l=1; l<=layers; l++) 
	icone_stencil.push(coupling_nodes_table.e(cn_indx2,l)+1);
    }
    // Add fictitious node to stencil
    r = coupl2fic.find(node);
    assert(r!=re);
    icone_stencil.push(r->second+1);
    nelem_nlr++;
  }

#if 0
  printf("nx: %d, nelem_nlr %d, size of icone %d\n",nx,nelem_nlr,icone_stencil.size());
  for (int j=0; j<nelem_nlr; j++) {
    printf("e=%d:",j);
    for (int l=0; l<nx+1; l++) 
      printf(" %d",icone_stencil.e(j,l));
    printf("\n");
  }
#endif

  null_vort *nv = new null_vort;
  Elemset *elemset = nv;
  elemset->type = local_copy("null_vort");
  elemset->nelem = nelem_nlr; 
  elemset->nel   = nx+1; 
  elemset->elem_conne = new int[nx+1];
  elemset->ndof  = dofmap->ndof; 
  elemset->nelprops = 0; 
  elemset->neliprops = 0; 
  elemset->nelprops_add = 0; 
  elemset->neliprops_add = 0; 
  elemset->elemprops_add = 0; 
  elemset->elemiprops_add = 0; 
  elemset->thash = thash; 
  elemset->icone = icone_stencil.buff(); 
  elemset->elemprops  = NULL; 
  elemset->elemiprops  = NULL; 
  elemset->elem_prop_names  = NULL; 
  elemset->epart = NULL;
  elemset->isfat = 0;

  // Registers the elemset in the static table.
  // This should be done for all the elemsets. 
  string name = Elemset::anon;
  ::get_string(thash,"name",name,1);
  elemset->register_name(name,elemset->type);
  thash->register_name(elemset->name());

  // Append to elemsetlist. This is the list of
  // elemsets to be processed each time we call `assemble'. 
  da_append(mesh->elemsetlist,&elemset);

  nv->fic_dof = fic_dof;

  icone.clear();
  coupling_nodes_table.clear();
  coupling_nodes_map.clear();
  graph.clear();
  coupl2fic.clear();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void null_vort::lag_mul_dof(int jr,int &node,int &dof) {
  assert(jr==1);
  node = nel;
  dof = fic_dof;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void null_vort::init() {
  int ierr;
  //o The number of spatial dimensions
  TGETOPTDEF_ND(thash,int,ndim,0);
  // Number of nodes in the stencil. (The -1 corresponds to
  // the fictitious node)
  nx = nel-1;
  // Build the cloud 
  // compute coefs. for d/dx and d/dy
  int derivs[] = {1,0,0,1};
  // use second order in both directions
  int npol[] = {2,2};
  // initialize the object
  cloud.init(ndim,nx,2,derivs,npol);
  xlocc.resize(2,nel,ndim);
  x0.resize(1,ndim);
  w.resize(2,nx,2);
  ww.resize(2,nel,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void null_vort::res(int k,FastMat2 &U,FastMat2 & r,FastMat2 & lambda,
		    FastMat2 & Jac) {
  // Make a local copy of node coordinates
  xlocc.set(xloc());
  // Coordinates of surface node. The first node in the stencil
  // is the center node (the node where we impose the null
  // vort. condition
  xlocc.ir(1,1);
  // Coordinates of all nodes
  x0.set(xlocc);
  xlocc.rs().is(1,1,nx);
  x.set(xlocc);
  xlocc.rs();
  x.rs();
  // Compute coefs.
  cloud.coef(x,w,x0);
  // warning if the computation was too bad conditioned
  double cond = cloud.cond();
  if (cond > 1e8) {
    PetscPrintf(PETSC_COMM_WORLD,
		"null_vort: object %d, bad conditioned cloud, cond %g\n",
		k,cond);
  }
  // Scale by volume of element
  w.scale(pow(cloud.h(),ndim));
  ww.set(0.).is(1,1,nx).set(w).rs();

  // Vorticity is w = dv/dx - du/dy
  // Coefs for dv/dx
  ww.ir(2,1);
  Jac.set(0.).ir(1,1).ir(3,2).set(ww);
  // Coefs for dv/dx
  ww.ir(2,2);
  Jac.ir(3,1).set(ww).scale(-1);
  ww.rs();
  Jac.rs();
#if 0
  xloc().print("");
  Jac.print("Jac:");
#endif

  // We compute the residual as Jac*U
  // r.prod(Jac,U,1,-1,-2,-1,-2).scale(-1); // should we scale by -1?
  r.prod(Jac,U,1,-1,-2,-1,-2);
  // We consider the Lag. Multiplier as symmetric, i.e. the
  // reaction is in the same direction as the restriction
  // so that we are `discarding' a linear combination of the
  // momentum equations for the nodes in the stencil
  lambda.set(0.).ctr(Jac,3,1,2);
}
