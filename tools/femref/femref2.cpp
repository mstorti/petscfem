//__INSERT_LICENSE__
// $Id: femref2.cpp,v 1.10 2004/12/14 14:38:44 mstorti Exp $

#include <string>
#include <list>
#include <iostream>
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
UniformMesh::
~UniformMesh() {
  for (int j=0; j<nelem; j++) {
    if (elem_ref.e(j)) delete elem_ref.e(j);
  }
  elem_ref.clear();
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
  elem_ref.resize(nelem);
  for (int j=0; j<nelem; j++) {
    elem_ref.e(j) = new ElemRef;
  }
}

// Simple macro-based STACK adaptor
// #define TOP(list_stack) (list_stack.begin())
// #define POP(list_stack) (list_stack.erase(list_stack.begin()))
// #define PUSH(list_stack,x) (list_stack.insert(list_stack.begin(),x))
// #define EMPTY(list_stack) (list_stack.size()==0)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define MAX_NODE 1000
void UniformMesh::
refine(RefineFunction f) {
  // The stack of geometrical objects while
  // traversing the tree
  list<GeomObject> go_stack;
  // The stack of splitters (nodes in the tree)
  list<ElemRef::iterator> split_stack;
  // The stack of indices on each splitting node
  list<int> split_indx_stack;
  // The stack of node indices 
  dvector<int> go_nodes;
  // The stack of positions in go_nodes
  list<int> go_nodes_pos;
  // Initialize first `refined' node
  last_ref_node = nnod;

  list<GeomObject>::iterator w,ws;
  // Loop over elements
  for (int k=0; k<nelem; k++) {
    // Clear all stacks
    go_stack.clear();
    split_stack.clear();
    split_indx_stack.clear();
    ElemRef &etree = *elem_ref.e(k);
    ElemRef::iterator q = etree.begin();

    go_stack.push_front(GeomObject());
    w = go_stack.begin();
    w->init(tmpl->type,&connec.e(k,0));
    w->make_canonical();
    w->print();
    split_stack.push_front(etree.begin());
    split_indx_stack.push_front(0);
    const int *w_nodes = w->nodes();
    for (int j=0; j<w->size(); j++)
      go_nodes.push(w_nodes[j]);
    go_nodes_pos.push_front(0);
    const Splitter *s = NULL;
    bool done = false;

    while(!done) {
      w = go_stack.begin();
      // here visit w...
      ElemRef::iterator q, qs, qs2, qfather;
      q = split_stack.front();
      int j = split_indx_stack.front();
      if (q != etree.end()) {
	s = q->splitter;
	qs = q.lchild();
	if (qs!=etree.end() && qs->so_indx==0) 
	  qs2 = qs;
	else qs2 = etree.end();
	split_stack.push_front(qs2);
	go_stack.push_front(GeomObject());
	ws = go_stack.begin();
	// Build `ws' from GeomObject `w' (parent) and splitter `s'
	// and subobject index `j'
	set(*w,s,j,*ws);
	ws->make_canonical();
	split_indx_stack.push_front(0);
      } else {
	// `q' is a leave for GO's (sure it isn't
	// a regular node for splitters)
	// Try to find a right sibling, or a father
	// that has a right sibling
	while (true) {
	  // Check if we are at the root
	  if (split_stack.size()<=1) {
	    done = true; break;
	  }
	  qfather = q;
	  qfather++;
	  assert(qfather != etree.end());
	  s = qfather->splitter;

	  go_stack.pop_front();
	  split_stack.pop_front();
	  split_indx_stack.pop_front();
	  if (j<s->size()) {
	    int jsib = j+1;
	    go_stack.push_front(GeomObject());
	    ws = go_stack.begin();
	    // Build `ws' from GeomObject `w' (parent) and splitter `s'
	    // and subobject index `jsib'
	    set(*w,s,jsib,*ws);
	    ws->make_canonical();

	    // Find next node on the splitting tree or end()
	    qs = qfather.lchild();
	    while (qs != etree.end()) {
	      if (qs->so_indx == jsib) break;
	      qs++;
	    }

	    // Push new state in the stacks
	    split_stack.push_front(qs);
	    split_indx_stack.push_front(jsib);
	    break;
	  }
	}
      }
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
  while (q!=etree.end()) {
    // Each node of the `etree' is a splitter!!
    // create the nodes needed by this splitting
    const Splitter *split = *q;
    int sz = split->size();
    for (int k=0; k<sz; k++) {
      GeomObject::Type t;
      const int *so_local_nodes = split->nodes(k,t);
      GeomObject sgo(t);
      int sgo_sz = sgo.size();
      sgo_nodes.resize(sgo_sz);
      for (int k=0; k<sgo_sz; k++) {
	int n1 = so_local_nodes[2*k];
	if (n1<w->size()) 
	  sgo_nodes.e(k) = w_nodes[k];
	else {
	  n2 = so_local_nodes[2*k];
	  n = n1*MAX_NODE+n2;
	  map<int,int>::iterator q = hash2node.find(n);
	  int ref_node;
	  if (q==hash2node.end()) {
	    ref_node = last_ref_node++;
	    hash2node[n] = ref_node;
	  } else {
	    ref_node = q->second;
	  }
	  sgo_nodes.e(k) = ref_node;
	}
      }
    }
    q = q.lchild();
  }
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::
set(const GeomObject &go,const Splitter *s,
    int indx,GeomObject &sgo) {
  const int *go_nodes = go.nodes();
  GeomObject::Type t;
  const int *local_nodes = s->nodes(indx,t);
  GeomObject::Template *tmpl 
    = GeomObject::get_template(t);
  int 
    so_sz = tmpl->size_m,
    sz = go.size();
  dvector<int> sgo_nodes;
  sgo_nodes.mono(so_sz);
  for (int k=0; k<so_sz; k++) {
    int ln1 = local_nodes[2*k];
    int n1 = go_nodes[ln1];
    if (ln1>=sz) {
      int ln2 = local_nodes[2*k+1];
      int n2 = go_nodes[ln2];
      int node_hash = combine(n1,n2);
      map<int,int>::const_iterator it 
	= hash2node.find(node_hash);
      if (it != hash2node.end()) {
	n1 = it->second;
      } else {
	n1 = last_ref_node++;
	hash2node[node_hash] = n1;
      }
    }
    sgo_nodes.e(k) = n1;
  }
  sgo.init(t,sgo_nodes.buff());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int UniformMesh::
combine(int n1,int n2) const {
  if (n2<n1) return combine(n2,n1);
  else return n1*MAX_NODE+n2;
}
