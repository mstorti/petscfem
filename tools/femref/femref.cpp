//__INSERT_LICENSE__
// $Id: femref.cpp,v 1.7 2004/11/21 23:43:09 mstorti Exp $

#include <string>
#include <limits.h>

using namespace std;

#include "./femref.h"

GeomObject::GeomObject(Type t,const int *nodes_a) 
  :  canonical(0), go_template(NULL) {

#define SET_TEMPL(CLASS) \
else if(t==CLASS##T) go_template = &CLASS##Template

  if (t==NULL_TYPE) assert(0);
  // SET_TEMPL(OrientedEdge);
  SET_TEMPL(Edge);
  // SET_TEMPL(OrientedTri);
  SET_TEMPL(Tri);
  SET_TEMPL(OrientedTetra);
  else assert(0);

  int sz = size();
  nodes.mono(sz);
  cs = 0;
  for (int j=0; j<sz; j++) {
    int node = nodes_a[j];
    cs += node;
    nodes.ref(j) = node;
  }

}

void GeomObject::make_canonical() {
  // iterate through all possible perms. and
  // remain with lowest (in lexicographical order)
  // numbering
  int np = nperms();
  int sz = size();
  dvector<int> aux;
  aux.mono(sz);
  int *cmin = &aux.ref(0);
  int *nodesp = &nodes.ref(0);
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

bool GeomObject::equal(GeomObject &go) {
  if (type() != go.type()) return false;
  if (csum() != go.csum()) return false;
  make_canonical();
  go.make_canonical();
  
}

int GeomObject::NULL_NODE = INT_MAX;

GeomObject::Template
::Template(int sz,GeomObject::Type t,
	     int dim_a,int nperms_a,const int *perms) 
  : size_m(sz), type(t), dim_m(dim_a), nperms_m(nperms_a) {
  perms_v.resize(nperms_m*size_m);
  const int *q = perms;
  for (int j=0; j<nperms_m*size_m; j++)
    perms_v.e(j) = *q++;
  assert(*q==GeomObject::NULL_NODE);
  perms_v.reshape(2,nperms_m,size_m);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static int EdgeTemplatePerm_v[] = {1,0,GeomObject::NULL_NODE};
GeomObject::Template 
GeomObject::EdgeTemplate(2,GeomObject::EdgeT,
			 1,1,EdgeTemplatePerm_v);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static int TriTemplatePerm_v[] = {1,2,0,
				  2,0,1,
				  0,2,1,
				  2,1,0,
				  1,0,2,GeomObject::NULL_NODE};
GeomObject::Template 
GeomObject::TriTemplate(3,GeomObject::TriT,
			 2,5,TriTemplatePerm_v);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static int OrientedTetraTemplate_v[] = {1,2,0,3,
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
GeomObject::Template 
GeomObject::OrientedTetraTemplate(4,GeomObject::OrientedTetraT,
				  3,11,OrientedTetraTemplate_v);

int main() { 
  int v1[] = {0,1};
  GeomObject edge(GeomObject::EdgeT,v1);
  int v2[] = {23,35,2};
  GeomObject tri(GeomObject::TriT,v2);
  tri.make_canonical();
  dvector<int> nodes;
  nodes.mono(4);
  nodes.ref(0) = 34;
  nodes.ref(1) = 98;
  nodes.ref(2) = 91;
  nodes.ref(3) = 76;
  GeomObject tetra(GeomObject::OrientedTetraT,nodes.buff());
  tetra.make_canonical();
#if 0
  while(1) {
    for (int j=0; j<4; j++) 
      nodes.ref(j) = rand() % 100;
    GeomObject tetra(GeomObject::OrientedTetraT,nodes.buff());
    tetra.make_canonical();
  }
#endif
}
