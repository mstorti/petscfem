//__INSERT_LICENSE__
// $Id: femref.cpp,v 1.24 2004/12/19 16:17:44 mstorti Exp $

#include <string>
#include <list>
#include <limits.h>
#include "./hasher.h"

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main1() { 
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
  list<Mesh::iterator> its;
  for (int j=0; j<5; j++) {
    for (int k=0; k<4; k++) {
      it.set(j,GeomObject::OrientedTriT,k);
      mesh.set(it,go);
      go.print();
#if 0
      it = mesh.find(go);
      if (mesh.is_end(it)) printf("not found\n");
      else printf("found at elem %d, position %d\n",it.obj,it.subobj);
#else
      mesh.find(go,its);
      if (its.begin()==its.end()) 
	printf("not found\n");
      else {
	printf("found at: ");
	list<Mesh::iterator>::iterator q = its.begin();
	while (q!=its.end()) {
	  printf("(elem %d, position %d) ",q->obj,q->subobj);
	  q++;
	}
	printf("\n");
      }
#endif
    }
  }
}

double
rf(GeomObject &go,const double *xnod) {
  return 0.1;
}

int main() { 

  Tetra2TetraSplitterClass 
    Tetra2TetraSplitter2;

  UniformMesh mesh(OrientedTetraTemplate,3);
  mesh.read("tetra.nod","tetra.con");
  // mesh.refine(rf);
  UniformMesh::visitor vis;
  vis.init(mesh,0);
  while (!vis.so_end()) {
    GeomObject &go = vis.ref_stack.front().go;
    // go.print();
    if (vis.is_leave() && vis.ref_level()<=1)
      vis.refine();
    vis.so_next();
  }
}
