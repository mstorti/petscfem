//__INSERT_LICENSE__
// $Id: femref.cpp,v 1.30 2004/12/19 23:16:04 mstorti Exp $

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
class NodeInfoSum : public NodeInfo {
public:
  int sum;
  ~NodeInfoSum() { }
  void print() { 
    printf("NodeInfo %p sum %d\n",
	   this,sum); 
  }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
NodeInfo*
node_sum(const NodeInfo &ni1,
	 const NodeInfo &ni2) {
  const NodeInfoSum *ni1_p 
    = dynamic_cast<const NodeInfoSum *>(&ni1);
  const NodeInfoSum *ni2_p 
    = dynamic_cast<const NodeInfoSum *>(&ni2);

  NodeInfoSum *ni12_p = new NodeInfoSum;
  ni12_p->sum = ni1_p->sum + ni2_p->sum;
  return ni12_p;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int main() { 

  UniformMesh mesh(OrientedTetraTemplate,3);
  mesh.read("tetra.nod","tetra.con");
  // mesh.refine(rf);
  UniformMesh::visitor vis;

  // Refine
  vis.init(mesh);
  while (!vis.end()) {
    if (vis.is_leave() && vis.ref_level()<=0)
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
  vis.node_comb_fun = &node_sum;
  vis.init(mesh);
  while (vis.next(1)) {  }
}
