//__INSERT_LICENSE__
// $Id: femref2.cpp,v 1.1 2004/12/05 15:38:37 mstorti Exp $

#include <string>
#include <list>
#include <limits.h>
#include "./hasher.h"

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

MD5Hasher hasher;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int OrientedTetraTemplateClass
::perm_v[] = {1,2,0,3,
	      2,0,1,3,
	      0,3,1,2,
	      1,0,3,2,
	      3,1,0,2,
	      1,3,2,0,
	      3,2,1,0,
	      2,1,3,0,
	      0,2,3,1,
	      2,3,0,1,
	      3,0,2,1,GeomObject::NULL_NODE};

int 
OrientedTetraTemplateClass 
::faces[] = {0,1,3,
	     1,2,3,
	     2,0,3,
	     0,2,1,GeomObject::NULL_NODE};

OrientedTetraTemplateClass
OrientedTetraTemplate;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int OrientedTriTemplateClass
::perm_v[] = {1,2,0,
	      2,0,1,
	      0,2,1,
	      2,1,0,
	      1,0,2,GeomObject::NULL_NODE};

OrientedTriTemplateClass
OrientedTriTemplate;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Tetra2TetraSplitterClass 
Tetra2TetraSplitter;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GeomObject::GeomObject(Type t,const int *nodes_a) 
  :  canonical(0), go_template(NULL) { init(t,nodes_a); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void GeomObject::init(Type t,const int *nodes_a) {

  clear();

  // Set template ptr from type
  go_template = get_template(t);
  // resize node array
  int sz = size();
  nodes_m.mono(sz);
  // Copy node array in internal vector
  // and compute check-sum
  if (nodes_a) {
    hasher.reset();
    for (int j=0; j<sz; j++) {
      int node = nodes_a[j];
      hasher.hash(node);
      nodes_m.ref(j) = node;
    }
    cs = hasher.val();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GeomObject::Template* GeomObject::get_template(Type t) {

#define SWITCH(TYPE) \
 case TYPE##T: return &TYPE##Template

  switch(t) {
  case NULL_TYPE: return NULL;
    // SWITCH(OrientedEdge);
    // SWITCH(Edge);
    SWITCH(OrientedTri);
    // SWITCH(Tri);
    SWITCH(OrientedTetra);
  default: assert(0);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void GeomObject::make_canonical() {
  // Iterate through all possible perms. and
  // remain with lowest (in lexicographical order)
  // numbering
  
  // Apply all permutations and compare
  // with current minimum.
  // The comparison is done at the same time
  // the permutation is being built, so that
  // if the first node is higher than the
  // first node of th current minimum then
  // we skip to the next permutation. 

  // Number of permutations for current geometry
  int np = nperms();
  // Number of nodes in object
  int sz = size();
  // aux/cmin will be the current minimum
  // ordering of nodes
  dvector<int> aux;
  aux.mono(sz);
  int *cmin = &aux.ref(0);
  // Ptr to internal nodes
  int *nodesp = &nodes_m.ref(0);
  //#define DBG
#ifdef DBG
  printf("entered: ");
  for (int kk=0; kk<sz; kk++) 
    printf(" %d",nodesp[kk]);
  printf("\n");
#endif
  // Make current minimum = current ordering
  for (int kk=0; kk<sz; kk++) 
    cmin[kk] = nodesp[kk];
  // Loop over permutations that leave
  // invariant this shape
  for (int j=0; j<np; j++) {
    // Go constructing the permutation. 
    // If it is higher skip. Else replace current min
    const int *permv = perm(j);
#if 0 && defined DBG
    printf("perm %d: ",j);
#endif
    for (int k=0; k<sz; k++) {
      int perm_node = nodesp[permv[k]];
      int node = cmin[k];
      if (perm_node < node) {
	// Permuted numeration is lower (lexicographically)
	for (int kk=0; kk<sz; kk++) 
	  cmin[kk] = nodesp[permv[kk]];
	break;
      } else if (perm_node>node) 
	// Permuted numeration is higher (lexicographically)
	// -> continue with next permutation
	break;
    }
#if 0 && defined DBG
    for (int kk=0; kk<sz; kk++) 
      printf(" %d",cmin[kk]);
    printf("\n");
#endif
  }
  // Copy lowest ordering to internal ordering
  for (int kk=0; kk<sz; kk++) 
    nodesp[kk] = cmin[kk];
#if 0
  printf("canonical form: ");
  for (int kk=0; kk<sz; kk++) 
    printf(" %d",nodesp[kk]);
  printf("\n");
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool GeomObject::equal(GeomObject &go) {
  // Objects should have the same type and check-sum
  if (type() != go.type()) return false;
  if (csum() != go.csum()) return false;
  // Bring both of them to canonical ordering
  make_canonical();
  go.make_canonical();
  // Compare node sequences
  int *np = nodes_m.buff();
  int *gonp = go.nodes_m.buff();
  for (int j=0; j<go_template->size_m; j++)
    if (*np++ != *gonp++) return false;
  return true;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int GeomObject::NULL_NODE = INT_MAX;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GeomObject::Template
::Template(int sz,GeomObject::Type t,
	   int dim_a,int nperms_a,const int *perms,
	   const char *label_a) 
  : size_m(sz), type(t), dim_m(dim_a), nperms_m(nperms_a),
    label(label_a) {
  // Copy external array to internal
  perms_v.resize(nperms_m*size_m);
  const int *q = perms;
  for (int j=0; j<nperms_m*size_m; j++)
    perms_v.e(j) = *q++;
  // Check end of sequence
  assert(*q==GeomObject::NULL_NODE);
  // reshape internal array
  perms_v.reshape(2,nperms_m,size_m);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void GeomObject::print(const char*s) const {
  if (!go_template) {
    printf("<empty>\n");
    return;
  }
  if(s) printf("%s: ",s);
  printf("<%s> ",go_template->label);
  for (int j=0; j<go_template->size_m; j++) 
    printf("%d ",nodes_m.ref(j));
  printf("\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void GeomObject
::set(Type t,int j,GeomObject &go) const {
  // Set #go# to internal subobject #j of type #t#. 
  // Template of subobject
  Template * go_tmpl = get_template(t);
  // Size (number of nodes) of subobject
  int go_sz = go_tmpl->size_m;
  // List of nodes for internal subobject
  dvector<int> go_nodes;
  go_nodes.mono(go_sz);
  // Nodes of subobject (local to large object)
  const int *local_nodes = go_template->nodes(t,j);
  const int *my_nodes = nodes();
  // Define sub-object nodes
  for (int k=0; k<go_sz; k++) 
    go_nodes.e(k) = my_nodes[local_nodes[k]];
  // Create sub-object
  go.init(t,go_nodes.buff());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int GeomObject::find(GeomObject &go) const {
  // Find sub-object go in large object *this. Returns
  // number of subobject in large object. Pay attention to
  // the fact that the sub-object numbering is relative to
  // the current large-object node numbering. 

  // Shape of sub-object
  Type typ = go.type();
  // Nbr of nodes in subobject
  int n_so = size(typ);
  // Internal sub-object of the given type
  GeomObject subobj;
  for (int j=0; j<n_so; j++) {
    // create j-th sub-object
    set(typ,j,subobj);
    // compare objects
    if (go.equal(subobj)) return j;
  }
  return -1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::set(iterator it,GeomObject &go) {
  // assert(it.t == tmpl->type);
  // check that large-object number is OK.
  assert(0 <= it.obj && it.obj<nelem);
  // Number of sub-objects of type t in large-object
  int sz = tmpl->size(it.t);
  assert(it.subobj < sz);
  // Shape of the sub-object
  GeomObject::Template *so_tmpl = 
    GeomObject::get_template(it.t);
  // Nodes in local numbering
  const int *local_nodes = tmpl->nodes(it.t,it.subobj);
  // Nodes of the subobject
  dvector<int> so_nodes;
  int so_sz = so_tmpl->size_m;
  so_nodes.mono(so_sz);
  for (int j=0; j<so_sz; j++)
    so_nodes.e(j) = connec.e(it.obj,local_nodes[j]);
  go.init(it.t,so_nodes.buff());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void 
UniformMesh::find(GeomObject &go,list<iterator> &its) {
  // This is a wrapper, returns the list
  iterator it;
  its.clear();
  find(go,&its,it);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Mesh::iterator 
UniformMesh::find(GeomObject &go) { 
  // This is a wrapper, returns the first iterator
  iterator it;
  find(go,NULL,it);
  return it;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::
find(GeomObject &go,list<iterator> *its,iterator &it) {
  // Finds the elements connected to #go#
  // by performing the intersection of the sets of
  // large-objects connected to the nodes in #go#. 
  int sz = go.size();
  assert(sz>0);
  // `p' is an array of cursors. `[p(j),pe(j))' is the 
  // set of elements connected to node `j' in `go'. 
  dvector<int> p, pe;
  p.mono(sz); pe.mono(sz);
  const int *nodes_p = go.nodes();
  int emin, emax, jmin;
  go.print("searching face: ");
  for (int j=0; j<sz; j++) {
    int node = nodes_p[j];
    p.e(j) = n2e_ptr.e(node);
    pe.e(j) = n2e_ptr.e(node+1);
    // Normally, each node should have at least one element
    // connected to it
    assert(pe.e(j)>p.e(j));
  }
#if 0
  for (int j=0; j<sz; j++) {
    printf("ngbrs elems for node %d: ",nodes_p[j]);
    for (int k=p.e(j); k<pe.e(j); k++) {
      printf("%d ",n2e.e(k));
    }
    printf("\n");
  }
#endif
  while (1) {
    // Search for min and max of elements pointed
    // by current position of cursors. 
    // printf("elems in ptrs: ");
    for (int j=0; j<sz; j++) {
      int pp = p.e(j);
      int ele = n2e.e(pp);
      // printf("%d ",ele);
      if (j==0) emin = emax = ele; 
      if (ele < emin) emin = ele;
      else if (ele > emax) emax = ele;
    }
    // printf("\n");
    // If all elems are equal, then this element
    // is a candidate to contain the object
    int indx;
    if (emin==emax) {
      // printf("candidate elem %d\n",emin);
      // Construct the large-object
      GeomObject super(tmpl->type,&connec.e(emin,0));
      // Look for `go' in this large-object
      indx = super.find(go);
      if (indx>=0) {
	// If `its==NULL' then the value is returned via `it',
	// else it is appended to `*its' and process continue. 
	if (!its) {
	  it = iterator(emin,go.type(),indx);
	  return;
	} else {
	  its->push_back(iterator(emin,go.type(),indx));
	}
      }
    } 
    // Advance all cursors that are pointing to the lowest
    // element.
    for (int j=0; j<sz; j++) {
      int pp = p.e(j);
      if (n2e.e(pp) == emin) {
	// We reached the end for at least one
	// of the list of elements. 
	if (++p.e(j) >= pe.e(j)) {
	  // If we are returning the list of iterators,
	  // then return the empty list. If not, return
	  // the empty iterator. 
	  if (!its) it = iterator();
	  return;
	} 
      }
    }
  }
  // we should never reach here
  assert(0); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Generates a random choice of `N' elements
// from the range `[0,M)', with `M<=N'. If `M=N'
// then this is a random permutation of the integers in `[0,N)'. 
void rand_perm(vector<int> &perm,int N,int M=-1) {
  if (M<0) M=N;
  vector<int> aux;
  aux.resize(M);
  perm.resize(N);
  for (int j=0; j<M; j++) aux[j] = j;
  for (int j=0; j<N; j++) {
    int k = rand() % N;
    perm[j] = aux[k];
    for (int kk=k+1; kk<M-j; kk++) 
      aux[kk-1] = aux[kk];
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Ctor from dimensions and shape
UniformMesh::
UniformMesh(GeomObject::Template &tmpl_a,int ndim_a) 
  : tmpl(&tmpl_a), nel(tmpl_a.size_m), ndim(ndim_a) { 
  // reshape dvectors to specicied shape
  coords.reshape(2,0,ndim);
  connec.reshape(2,0,nel);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Read the mesh from specified files
void UniformMesh::
read(const char *node_file,
     const char *conn_file) {
  // This is used auxiliary
  LinkGraph lgraph;
  // Reads coors file and resize
  coords.cat(node_file).defrag();
  // Reads connectivity file and resizes
  connec.cat(conn_file).defrag();
  // Gets number of nodes
  nnod = coords.size(0);
  // Gets number of elements per node.
  // This comes from the specified shape. 
  nel = connec.size(1);
  // Add connections to graph
  lgraph.init(nnod);
  nelem = connec.size(0);
  for (int ele=0; ele<nelem; ele++) {
    for (int k=0; k<nel; k++) {
      int node = connec.e(ele,k);
      // printf("add node %d, elem %d\n",node,ele);
      lgraph.add(node,ele);
    }
  }

  GSet ngbrs;
#if 0
  for (int node=0; node<nnod; node++) {
    // printf("node %d, elems ",node);
    // For each node get list of ngbrs and
    ngbrs.clear();
    lgraph.set_ngbrs(node,ngbrs);
    GSet::iterator p = ngbrs.begin();
    while (p!=ngbrs.end()) printf("%d ",*p++);
    printf("\n");
  }
#endif

  // Pass the connectivity in graph `lgraph' to
  // to (n2e,n2e_ptr)
  // Dimension `n2e_ptr'
  n2e_ptr.resize(nnod+1);
  // Initialize `n2e_ptr'
  for (int j=0; j<nnod; j++) n2e_ptr.ref(j) = 0;
  // Cummulate size of `ngbrs' in order
  // to define `n2e_ptr'
  int nadj = 0;
  for (int nod=0; nod<nnod; nod++) {
    n2e_ptr.e(nod) = nadj;
    ngbrs.clear();
    lgraph.set_ngbrs(nod,ngbrs);
    nadj += ngbrs.size();
    GSet::iterator p = ngbrs.begin();
    while (p!=ngbrs.end()) n2e.push(*p++);
  }
  // Set last (nnod+1) ptr
  n2e_ptr.ref(nnod) = nadj;
  // Free memory in `lgraph'
  lgraph.clear();

  const Splitter *p = NULL;
  elem_ref.resize(nelem);
}

double
rf(GeomObject &go,const double *xnod) {
  return 0.1;
}

void UniformMesh::
refine(RefineFunction f) {
  list<GeomObject> go_stack;
  for (int k=0; k<nelem; k++) {
    go_stack.clear();
  }
}
