//__INSERT_LICENSE__
// $Id: getsurf.cpp,v 1.9 2005/01/10 20:13:35 mstorti Exp $

#include <string>
#include <list>
#include <ctime>
#include <unistd.h>
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

int main(int argc,char **argv) {

  const char *iconef = NULL, *xnodf = NULL, 
    *statef = NULL, *sconf = NULL, *graduf = NULL;
  int base=0, nread; 
  char c;
  while ((c = getopt (argc, argv, "b:c:x:u:s:g:")) != -1) {
    switch (c) {
    case 'b':
      scanf(optarg,"%d",&base);
      assert(nread==1);
      break;
    case 'c':
      iconef = optarg;
      break;
    case 'x':
      xnodf = optarg;
      break;
    case 'u':
      statef = optarg;
      break;
    case 's':
      sconf = optarg;
      break;
    case 'g':
      graduf = optarg;
      break;
    default:
      abort ();
    }
  }

  assert(xnodf);
  assert(iconef);
  int ndof = 4, ndim = 3;
  time_t start, end;
  dvector<double> u,x;
  const GeomObject::Template *mesh_tmpl = &OrientedTetraTemplate;
  UniformMesh mesh(*mesh_tmpl,3);
  int mesh_nel = mesh.tmplt()->size_m;
  mesh.read(xnodf,iconef,1);
  u.reshape(2,0,ndof);
  u.cat(statef);
  u.defrag();

  x.reshape(2,0,ndim);
  x.cat(xnodf);
  x.defrag();
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
      int hash_val = inv_face.csum();
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
	int hv = face.csum();
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
  FastMat2 U(2,mesh_nel,ndof), 
    A(2,ndim+1,ndof), grad_U(2,ndim,ndof),
    X(2,ndim+1,mesh_nel), invX;
  X.ir(2,ndim+1).set(1.0).rs();
  FILE *fid = fopen(sconf,"w");
  FILE *fidgu = fopen(graduf,"w");
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
    if (VERBOSE) {
      printf("elem %d, ");
      go.print();
      printf("face %d, ",q->second.elem,q->second.face);
      face.print();
      printf("opposing node %d\n",opp_node);
    }
    X.is(2,1,ndim);
    for (int j=0; j<mesh_nel; j++) {
      int node;
      if (j<face_nel) node = face_nodes[j];
      else node = opp_node;
      U.ir(1,j+1).set(&u.e(node,0));
      X.ir(1,j+1).set(&x.e(node,0));
    }
    U.rs(); X.rs();
    invX.inv(X);
    A.prod(invX,U,1,-1,-1,2);
    A.is(1,1,ndim);
    grad_U.set(A);
    A.rs();
    // grad_U.print("");
    for (int j=0; j<face_nel; j++) {
      int node = face_nodes[j];
      fprintf(fid,"%d ",node);
    }
    fprintf(fid,"\n");
    double *buff = grad_U.storage_begin();
    for (int j=0; j<ndim*ndof; j++) 
      fprintf(fidgu,"%lg ",buff[j]);
    fprintf(fidgu,"\n");
  }
  fclose(fid);
  fclose(fidgu);
}
