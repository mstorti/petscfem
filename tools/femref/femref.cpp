//__INSERT_LICENSE__
// $Id: femref.cpp,v 1.15 2004/11/23 12:37:30 mstorti Exp $

#include <string>
#include <limits.h>

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

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
GeomObject::GeomObject(Type t,const int *nodes_a) 
  :  canonical(0), go_template(NULL) { init(t,nodes_a); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void GeomObject::init(Type t,const int *nodes_a) {

  clear();

#if 0
#define SET_TEMPL(CLASS) \
else if(t==CLASS##T) go_template = &CLASS##Template

  if (t==NULL_TYPE) assert(0);
  // SET_TEMPL(OrientedEdge);
  // SET_TEMPL(Edge);
  SET_TEMPL(OrientedTri);
  // SET_TEMPL(Tri);
  SET_TEMPL(OrientedTetra);
  else assert(0);
#endif

  go_template = get_template(t);
  int sz = size();
  nodes_m.mono(sz);
  if (nodes_a) {
    cs = 0;
    for (int j=0; j<sz; j++) {
      int node = nodes_a[j];
      cs += node;
      nodes_m.ref(j) = node;
    }
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
  // iterate through all possible perms. and
  // remain with lowest (in lexicographical order)
  // numbering
  int np = nperms();
  int sz = size();
  dvector<int> aux;
  aux.mono(sz);
  int *cmin = &aux.ref(0);
  int *nodesp = &nodes_m.ref(0);
#define DBG
#ifdef DBG
  printf("entered: ");
  for (int kk=0; kk<sz; kk++) 
    printf(" %d",nodesp[kk]);
  printf("\n");
#endif
  for (int kk=0; kk<sz; kk++) 
    cmin[kk] = nodesp[kk];
  for (int j=0; j<np; j++) {
    // Go constructing the permutation. 
    // If it is higher skip. Else replain current min
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
  for (int kk=0; kk<sz; kk++) 
    nodesp[kk] = cmin[kk];
#if 1
  printf("canonical form: ");
  for (int kk=0; kk<sz; kk++) 
    printf(" %d",nodesp[kk]);
  printf("\n");
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bool GeomObject::equal(GeomObject &go) {
  if (type() != go.type()) return false;
  if (csum() != go.csum()) return false;
  make_canonical();
  go.make_canonical();
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
  perms_v.resize(nperms_m*size_m);
  const int *q = perms;
  for (int j=0; j<nperms_m*size_m; j++)
    perms_v.e(j) = *q++;
  assert(*q==GeomObject::NULL_NODE);
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
  Template * go_tmpl = get_template(t);
  int go_sz = go_tmpl->size_m;
  dvector<int> go_nodes;
  go_nodes.mono(go_sz);
  const int *local_nodes = go_template->nodes(t,j);
  const int *my_nodes = nodes();
  for (int k=0; k<go_sz; k++) 
    go_nodes.e(k) = my_nodes[local_nodes[k]];
  go.init(t,go_nodes.buff());
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int GeomObject::find(GeomObject &go) const {
  Type typ = go.type();
  // Template *templ = get_template(t);
  int n_so = size(typ);
  GeomObject subobj;
  for (int j=0; j<n_so; j++) {
    set(typ,j,subobj);
    if (go.equal(subobj)) return j;
  }
  return -1;
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define TEMPLATE(TYPE,size,dim,nperms)				\
GeomObject::Template						\
GeomObject::TYPE##Template(size,GeomObject::TYPE##T,		\
			 dim,nperms,TYPE##TemplatePerm_v,#TYPE)

static int EdgeTemplatePerm_v[] = {1,0,GeomObject::NULL_NODE};
TEMPLATE(Edge,2,1,1);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static int 
TriTemplatePerm_v[] = {1,2,0,
		       2,0,1,
		       0,2,1,
		       2,1,0,
		       1,0,2,GeomObject::NULL_NODE};
TEMPLATE(Tri,3,2,5);
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UniformMesh::set(iterator it,GeomObject &go) {
  // assert(it.t == tmpl->type);
  assert(0 <= it.obj && it.obj<nelem);
  int sz = tmpl->size(it.t);
  assert(it.subobj < sz);
  // go.init(it.t);
  GeomObject::Template *so_tmpl = 
    GeomObject::get_template(it.t);
  const int *local_nodes = tmpl->nodes(it.t,it.subobj);
  // Nodes of the subobject
  dvector<int> so_nodes;
  int so_sz = so_tmpl->size_m;
  so_nodes.mono(so_sz);
  for (int j=0; j<so_sz; j++)
    so_nodes.e(j) = connec.e(it.obj,local_nodes[j]);
  go.init(it.t,so_nodes.buff());
}

Mesh::iterator 
UniformMesh::find(GeomObject &go) { 
  int sz = go.size();
  assert(sz>0);
  dvector<int> p, pe, vals;
  p.mono(sz); pe.mono(sz); vals.mono(sz);
  const int *nodes_p = go.nodes();
  int emin, emax, jmin;
  go.print("searching face: ");
  for (int j=0; j<sz; j++) {
    int node = nodes_p[j];
    p.e(j) = n2e_ptr.e(node);
    pe.e(j) = n2e_ptr.e(node+1);
    assert(pe.e(j)>p.e(j));
  }
  bool found = false, done = false;
#if 1
  for (int j=0; j<sz; j++) {
    printf("ngbrs elems for node %d: ",nodes_p[j]);
    for (int k=p.e(j); k<pe.e(j); k++) {
      printf("%d ",n2e.e(k));
    }
    printf("\n");
  }
#endif
  while (!done) {
    // Search for min and max elements
    printf("elems in ptrs: ");
    for (int j=0; j<sz; j++) {
      int pp = p.e(j);
      int ele = n2e.e(pp);
      printf("%d ",ele);
      if (j==0) emin = emax = ele; 
      if (ele < emin) emin = ele;
      else if (ele > emax) emax = ele;
    }
    printf("\n");
    // If all elems are equal 
    int indx;
    if (emin==emax) {
      printf("candidate elem %d\n",emin);
      // GeomObject(Type t,const int *nodes_a=NULL);
      GeomObject super(tmpl->type,&connec.e(emin,0));
      indx = super.find(go);
      if (indx>=0) {
	found = true;
	return iterator(emin,go.type(),indx);
      }
    } 
    for (int j=0; j<sz; j++) {
      int pp = p.e(j);
      if (n2e.e(pp) == emin) {
	if (++p.e(j) >= pe.e(j)) {
	  done = true;
	  return iterator();
	} 
      }
    }
  }
  assert(0); // we should never reach here
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
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
int main() { 
#if 0
  int v1[] = {0,1};
  GeomObject edge(GeomObject::EdgeT,v1);
  int v2[] = {23,35,2};
  GeomObject tri(GeomObject::TriT,v2);
  tri.make_canonical();
  GeomObject tetra1, tetra2;
  vector<int> nodes(4);
  while(1) {
    rand_perm(nodes,4,6);
    tetra1.init(GeomObject::OrientedTetraT,&nodes[0]);
    tetra1.print("tetra1");

    rand_perm(nodes,4,6);
    tetra2.init(GeomObject::OrientedTetraT,&nodes[0]);
    tetra2.print("tetra2");

    printf("tetra1==tetra2? %s\n",(tetra1.equal(tetra2) ? "yes" : "no"));
  }
#endif
  UniformMesh mesh(OrientedTetraTemplate,3);
  mesh.read("tetra.nod","tetra.con");
  GeomObject go;
  Mesh::iterator it;
  for (int j=0; j<5; j++) {
    for (int k=0; k<4; k++) {
      it.set(j,GeomObject::OrientedTriT,k);
      mesh.set(it,go);
      go.print();
      it = mesh.find(go);
      if (mesh.is_end(it)) printf("not found\n");
      else printf("found at elem %d, position %d\n",it.obj,it.subobj);
    }
  }
}
