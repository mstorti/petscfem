//__INSERT_LICENSE__
// $Id: getsurf.cpp,v 1.5 2005/01/09 21:32:29 mstorti Exp $

#include <string>
#include <list>
#include <multimap.h>
// #include <algorithm>
#include <limits.h>
#include "./hasher.h"

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

struct FaceIterator {
  int elem;
  int face;
};

typedef multimap<int,FaceIterator> face_table_t;
typedef pair<int,FaceIterator> face_table_pair_t;

int main() { 

  UniformMesh mesh(OrientedTetraTemplate,3);
#if 0
#define DATA "/u/mstorti/PETSC/GARIBA/DATA"
  mesh.read(DATA "/proy.nod.tmp",DATA "/proy.con.tmp");
#else
  mesh.read("tetra.nod","tetra.con");
#endif
  UniformMesh::visitor vis;
  vis.visit_mode = UniformMesh::BreadthFirst;
  vis.init(mesh);
  face_table_t face_table;
  GeomObject face, inv_face;
  const GeomObject::Template *face_tmpl 
    = GeomObject::get_template(GeomObject::OrientedTriT);
  dvector<int> inds;
  int m = face_tmpl->size_m;
  inds.resize(m);
  int nint_faces=0;
  while (!vis.end()) {  
    GeomObject &go = vis.ref_stack.front().go;
    go.print();
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

      printf("face %d ",j);
      face.print();
      int hash_val = inv_face.hash_val();
      face_table_t::iterator q, 
	q1 = face_table.lower_bound(hash_val),
	q2 = face_table.upper_bound(hash_val);
      if (q1!=q2) {
	printf("possible collition with ");
	for (q = q1; q!=q2; q++)
	  printf("(elem %d, face %d) ", 
		 q->second.elem,q->second.face);
	printf("\n");
	nint_faces++;
      } else {
	FaceIterator fi;
	fi.elem = vis.elem_indx();
	fi.face = j;
	face_table
	  .insert(face_table_pair_t(face.hash_val(),fi));
      }
    }
    vis.next();
  }
  printf("%d internal faces, %d external faces\n",
	 nint_faces,face_table.size());
}
