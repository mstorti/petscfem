//__INSERT_LICENSE__
// $Id: getsurf.cpp,v 1.7 2005/01/09 23:53:58 mstorti Exp $

#include <string>
#include <list>
#include <multimap.h>
// #include <algorithm>
#include <limits.h>
#include "./hasher.h"

using namespace std;

#define VERBOSE 0

#include "./femref.h"
#include "./gtemplates.h"

struct FaceIterator {
  int elem;
  int face;
};

typedef multimap<int,FaceIterator> face_table_t;
typedef pair<int,FaceIterator> ft_pair_t;

int main() { 

  UniformMesh mesh(OrientedTetraTemplate,3);
#if 0
#define DATA "/u/mstorti/PETSC/GARIBA/DATA"
  mesh.read(DATA "/proy.nod.tmp",DATA "/proy.con.tmp");
#elif 0
  mesh.read("tetra.nod","tetra.con");
#else
  mesh.read("cube.nod.tmp","cube.con.tmp");
#endif
  UniformMesh::visitor vis, vis2;
  vis.visit_mode = UniformMesh::BreadthFirst;
  vis.init(mesh);
  vis2.init(mesh);
  face_table_t face_table;
  GeomObject face, inv_face, inv_face2;
  const GeomObject::Template *face_tmpl 
    = GeomObject::get_template(GeomObject::OrientedTriT);
  dvector<int> inds;
  int m = face_tmpl->size_m;
  inds.resize(m);
  int nint_faces=0;
  while (!vis.end()) {  
    GeomObject &go = vis.ref_stack.front().go;
    if (VERBOSE) {
      printf("elem %d, ",vis.elem_indx());
      go.print();
    }
    int nfaces = go.size(GeomObject::OrientedTriT);
    for (int j=0; j<nfaces; j++) {
      go.set(GeomObject::OrientedTriT,j,face);
      face.make_canonical();
      const int *nds = face.nodes();
      for (int l=0; l<m; l++) 
	inds.ref(l) = nds[l];
      // Invert face
      int x = inds.ref(1);
      inds.ref(1) = inds.ref(2);
      inds.ref(2) = x;
      inv_face.init(GeomObject::OrientedTriT,
		    inds.buff());
      inv_face.make_canonical();

      if (VERBOSE) {
	printf("face %d ",j);
	face.print();
      }
      int hash_val = inv_face.hash_val();
      face_table_t::iterator q, qq,
	q1 = face_table.lower_bound(hash_val),
	q2 = face_table.upper_bound(hash_val);
      bool add_face = true;
      if (q1!=q2) {
	if (VERBOSE) {
	  printf("possible collision with ");
	  for (q = q1; q!=q2; q++) 
	    printf("(elem %d, face %d) ", 
		   q->second.elem,q->second.face);
	  printf("\n");
	}
	int nfaces=0;
	for (q = q1; q!=q2; q++) {
	  vis2.init(q->second.elem);
	  GeomObject &go2 = vis2
	    .ref_stack.front().go;
	  go2.set(GeomObject::OrientedTriT,q->second.face,inv_face2);
	  if (inv_face.equal(inv_face2)) {
	    nfaces++;
	    qq = q;
	  }
	}
	if (nfaces!=1) {
	  printf("possible collision with ");
	  for (q = q1; q!=q2; q++) 
	    printf("(elem %d, face %d) ", 
		   q->second.elem,q->second.face);
	  printf("\n");
	}
	assert(nfaces==1 || nfaces==0);
	if (nfaces==1) {
	  face_table.erase(qq);
	  nint_faces++;
	  add_face = false;
	} 
      }
      if (add_face) {
	FaceIterator fi;
	fi.elem = vis.elem_indx();
	fi.face = j;
	int hv = face.hash_val();
	if (VERBOSE) {
	  printf("inserting face: hash %d, elem %d, j %d, ",
		 hv,fi.elem,fi.face);
	  face.print("");
	}
	face_table.insert(ft_pair_t(hv,fi));
      }
    }
    vis.next();
  }
  printf("%d total faces, %d internal, %d external\n",
	 nint_faces*2+face_table.size(),
	 nint_faces,face_table.size());
}
