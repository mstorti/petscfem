//__INSERT_LICENSE__
// $Id: getsurf.cpp,v 1.8 2005/01/10 16:12:18 mstorti Exp $

#include <string>
#include <list>
#include <ctime>
#include <multimap.h>
// #include <algorithm>
#include <limits.h>
#include "./hasher.h"
#include <src/fastmat2.h>

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

#define MAXH 10000

int main() { 

  time_t start, end;
  dvector<double> u;
  const GeomObject::Template *mesh_tmpl = &OrientedTetraTemplate;
  UniformMesh mesh(*mesh_tmpl,3);
  int mesh_nel = mesh.tmplt()->size_m;
#if 0
#define DATA "/u/mstorti/PETSC/GARIBA/DATA"
  mesh.read(DATA "/proy.nod.tmp",DATA "/proy.con.tmp");
#elif 0
  mesh.read("tetra.nod","tetra.con");
#else
  mesh.read("cube.nod.tmp","cube.con.tmp");
  u.cat("cube.state.tmp");
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
  int face_nel = face_tmpl->size_m;
  inds.resize(face_nel);
  int nint_faces=0;
  int tries=0;
  double table_size=0.0;
  int collis = 0;
  start = time(NULL);
  while (!vis.end()) {  
    GeomObject &go = vis.ref_stack.front().go;
    if (VERBOSE) {
      printf("elem %d, ",vis.elem_indx());
      go.print();
    }
    int nfaces = go.size(GeomObject::OrientedTriT);
    for (int j=0; j<nfaces; j++) {
      tries++;
      table_size += face_table.size();
      go.set(GeomObject::OrientedTriT,j,face);
      face.make_canonical();
      const int *nds = face.nodes();
      for (int l=0; l<face_nel; l++) 
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
      int hash_val = inv_face.csum() % MAXH;
      face_table_t::iterator q, qq,
	q1 = face_table.lower_bound(hash_val),
	q2 = face_table.upper_bound(hash_val);
      int nfaces=0, ncandi=0;
      for (q = q1; q!=q2; q++) {
	ncandi++;
	vis2.init(q->second.elem);
	GeomObject &go2 = vis2
	  .ref_stack.front().go;
	go2.set(GeomObject::OrientedTriT,q->second.face,inv_face2);
	if (inv_face.equal(inv_face2)) {
	  nfaces++;
	  qq = q;
	}
      }
      if (q1!=q2 && nfaces!=1 && VERBOSE) {
	printf("possible collision with ");
	for (q = q1; q!=q2; q++) 
	  printf("(elem %d, face %d) ", 
		 q->second.elem,q->second.face);
	printf("\n");
      }
      assert(nfaces<=1);
      if (nfaces==1) {
	// Face is duplicated (internal)
	face_table.erase(qq);
	nint_faces++;
      } else {
	// New face, add to table
	FaceIterator fi;
	fi.elem = vis.elem_indx();
	fi.face = j;
	int hv = face.csum() % MAXH;
	if (VERBOSE) {
	  printf("inserting face: hash %d, elem %d, j %d, ",
		 hv,fi.elem,fi.face);
	  face.print("");
	}
	face_table.insert(ft_pair_t(hv,fi));
      }
      collis += ncandi - nfaces;
    }
    vis.next();
  }
  printf("%d total faces, %d internal, %d external\n",
	 nint_faces*2+face_table.size(),
	 nint_faces,face_table.size());
  printf("tries %d, averg. table size %.2f, collisions %d\n",
	 tries,table_size/tries,collis);
  printf("face elimination %.2fsecs\n",
	 difftime(time(NULL),start));

  face_table_t::iterator q, qq,
    q1 = face_table.begin(),
    q2 = face_table.end();
  vector<int> node_mark;
  node_mark.resize(mesh_nel);
  for (q=q1; q!=q2; q++) {
    vis.init(q->second.elem);
    GeomObject &go = vis.ref_stack.front().go;
    go.set(GeomObject::OrientedTriT,q->second.face,face);
    for (int j=0; j<mesh_nel; j++) node_mark[j]=0;
    const int *face_nodes = face.nodes();
    const int *elem_nodes = go.nodes();
    for (int j=0; j<face_nel; j++) {
      for (int l=0; l<mesh_nel; l++) {
	if (elem_nodes[l] == face_nodes[j]) {
	  node_mark[l]=1; break;
	}
      }
    }
    int nopp=0, opp_node;
    for (int l=0; l<mesh_nel; l++) {
      if (node_mark[l]==0) {
	nopp++; opp_node = elem_nodes[l];
      }
    }
    printf("elem %d, ");
    go.print();
    printf("face %d, ",q->second.elem,q->second.face);
    face.print();
    printf("opposing node %d\n",opp_node);
  }
}
