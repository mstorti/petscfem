//__INSERT_LICENSE__
// $Id: nullvort.cpp,v 1.8 2003/02/28 15:37:38 mstorti Exp $

#include <src/nullvort.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/surf2vol.h>
#include <src/util2.h>
#include <src/linkgraph.h>
#include <src/cloud2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
null_vort::null_vort() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
null_vort::~null_vort() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void null_vort::read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap) {
  // Read options from data file
  thash.read(fstack);
  int ierr;

  // get options. Other options are processed in Surf2vol::factory
  //o The number of nodes per skin panel. 
  TGETOPTNDEF(&thash,int,nel_surf,0);
  //o The elemset that is from the fluid side. 
  TGETOPTDEF_S(&thash,string,volume_elemset,<none>);
  assert(volume_elemset!="<none>");

  // Call the `factory' (constructor) for the Surf2Vol object. 
  int identify_volume_elements, layers,
    use_exterior_normal,  ndimel;
  Surf2Vol *sv_gp_data=NULL;
  Elemset *vol_elem;
  Surf2Vol::factory(&thash, volume_elemset, 
		    nel_surf, sv_gp_data, 
		    vol_elem, identify_volume_elements, layers,
		    use_exterior_normal,  ndimel);
  int nel = nel_surf*(layers+1);

  dvector<int> icone;
  icone.reshape(2,0,nel);
  vector<string> tokens;
  vector<int> row;
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

  dvector<double> xnod;
  xnod.a_resize(2,mesh->nodedata->nnod,mesh->nodedata->nu);
  xnod.set(mesh->nodedata->nodedata);
  int ndim = mesh->nodedata->ndim;

  // Loop over nodes on the coupling surface
  GSet ngb;
  for (cn_indx=0; cn_indx<n_coupling_nodes; cn_indx++) {
    // Node number (0-based)
    int node = coupling_nodes_table.e(cn_indx,0);
    // get neighbors (on the surface)
    ngb.clear();
    graph.set_ngbrs(node,ngb);
    // Number of neighbors. We have a grid of
    // ngb.size() nodes on the surface. For each
    // node on the surface there is a row of (layers+1) nodes
    // (including the surface node) in the direction
    // nomal to the surface. This makes a total amount
    // of `n_cloud = ngb.size()*(layers+1)' nodes. 
    if (ngb.size()!=3) {
      printf("No 3 ngbrs node: %d\n",node);
      continue;
    }
    int nx = ngb.size() * (layers+1);

    // Build the cloud 
    Cloud2 cloud;
    int derivs[] = {1,0,0,1};
    int npol[] = {2,2};
    cloud.init(ndim,nx,2,derivs,npol);
    FastMat2 x(2,nx,ndim), x0(1,ndim), w(2,nx,2);
    // Coordinates of surface node
    x0.set(&xnod.e(node,0));
    // Coordinates of nodes in the cloud
    int k=0;
    GSet::iterator q, qe=ngb.end();
    for (q=ngb.begin(); q!=qe; q++) {
      int sf_node = *q;
      assert(coupling_nodes_map.find(sf_node)
	     !=coupling_nodes_map.end());
      cn_indx = coupling_nodes_map[sf_node];
      for (int j=0; j<=layers; j++) {
	// node in the cloud (0-based) (??)
	int node2 = coupling_nodes_table.e(cn_indx,j);
	x.ir(1,++k).set(&xnod.e(node2,0));
      }
    }
    x.rs();
    cloud.coef(x,w,x0);
    x0.print("x0: ");
    x.print("x: ");
    w.print("w: ");
  }

  icone.clear();
  xnod.clear();
  coupling_nodes_table.clear();
  coupling_nodes_map.clear();
  graph.clear();
}
