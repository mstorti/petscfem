//__INSERT_LICENSE__
// $Id: femref.cpp,v 1.36 2005/01/05 12:21:53 mstorti Exp $

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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double
rf(GeomObject &go,const double *xnod) {
  return 0.1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main() { 

  printf("starting main\n");
  UniformMesh mesh(OrientedTetraTemplate,3);
  mesh.read("tetra.nod","tetra.con");
  // mesh.refine(rf);
  UniformMesh::visitor vis;
  linear_combiner.mesh = &mesh;

  // Refine
  vis.init(mesh);
  while (!vis.end()) {
    if (vis.is_leave() && vis.ref_level()<=4)
      vis.refine(&Tetra2TetraSplitter);
    vis.next();
  }

#if 0
  // Print mesh
  vis.trace = 1;
  vis.init(mesh);
  while (!vis.end()) vis.next();

  // Print mesh down to level 1
  vis.trace = 1;
  vis.init(mesh);
  bool done;
  while (true) {
    if (vis.ref_level()<1) done = !vis.next();
    else done = !vis.level_next();
    if (done) break;
  }
#endif

  // Print mesh down to level 0
  vis.trace = 1;
  vis.node_comb = &linear_combiner;
  while (1) {
    vis.init(mesh);
    while (vis.next(3)) {  }
  }
}
