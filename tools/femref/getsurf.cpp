//__INSERT_LICENSE__
// $Id: getsurf.cpp,v 1.25 2005/01/16 19:59:02 mstorti Exp $

#include <string>
#include <list>
#include <set>
#include <ctime>
#include <unistd.h>
#include <multimap.h>
// #include <algorithm>
#include <limits.h>
#include "./hasher.h"
#include <src/fastmat2.h>

using namespace std;

#include "./femref.h"
#include "./gtemplates.h"

struct FaceIterator {
  int elem;
  int face;
};

typedef multimap<int,FaceIterator> face_table_t;
typedef pair<int,FaceIterator> ft_pair_t;

#define TRACE(j) printf("trace %d\n",j)

void getsurf(GetSurfCtx &ctx,
	     const dvector<int> &icone,
	     dvector<int> &surf_con,
	     dvector<int> &surf_nodes, 
	     int base, int verbose) {

  int nread; 
  char c;

  int ndof = 4, ndim = 3;
  ctx.ndim = ndim;
  time_t start, end;
  assert(!ctx.mesh);
  const GeomObject::Template *mesh_tmpl = &OrientedTetraTemplate;
  ctx.mesh = new UniformMesh(*mesh_tmpl,3);
  UniformMesh &mesh = *ctx.mesh;
  int mesh_nel = mesh.tmplt()->size_m;
  mesh.set_conn(icone,base);
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
    if (verbose) {
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

      if (verbose) {
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
      if (q1!=q2 && nfaces!=1 && verbose) {
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
	if (verbose) {
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
#if 0
  printf("%d total faces, %d internal, %d external\n",
	 nint_faces*2+face_table.size(),
	 nint_faces,face_table.size());
  printf("tries %d, averg. table size %.2f, collisions %d\n",
	 tries,table_size/tries,collis);
  printf("face elimination %.2fsecs\n",
	 difftime(time(NULL),start));
#endif

  face_table_t::iterator q, qq,
    q1 = face_table.begin(),
    q2 = face_table.end();
  vector<int> node_mark;
  node_mark.resize(mesh_nel);
  map<int,int> surf_nodes_map;
  int nfaces = face_table.size(),
    nsurf_nodes = 0;
  surf_con.a_resize(2,nfaces,face_nel);
  int surf_elem=0;

  for (q=q1; q!=q2; q++) {
    vis.init(q->second.elem);
    GeomObject &go = vis.ref_stack.front().go;
    go.set(GeomObject::OrientedTriT,q->second.face,face);
    for (int j=0; j<mesh_nel; j++) node_mark[j]=0;
    const int *face_nodes = face.nodes();
    const int *elem_nodes = go.nodes();
    for (int j=0; j<face_nel; j++) {
      int node = face_nodes[j];
      if (surf_nodes_map.find(node)
	  == surf_nodes_map.end()) {
	surf_nodes_map[node] = nsurf_nodes++;
	surf_nodes.push(node);
      }
      surf_con.e(surf_elem,j) 
	= surf_nodes_map[node];
      for (int l=0; l<mesh_nel; l++) {
	if (elem_nodes[l] == node) {
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
    assert(nopp=1);
    if (verbose) {
      printf("elem %d, ");
      go.print();
      printf("face %d, ",q->second.elem,q->second.face);
      face.print();
      printf("opposing node %d\n",opp_node);
    }
    surf_elem++;
  }
}

void comp_matrices(GetSurfCtx &ctx,
		   const dvector<int> &surf_con,
		   const dvector<int> &surf_nodes, 
		   dvector<double> &x, 
		   dvector<double> &surf_mass,
		   dvector<double> &node_mass,
		   int verbose) {

  printf("surf_con shape %d, %d\n",
	 surf_con.size(0),surf_con.size(1));
  printf("surf_nodes size %d\n",surf_nodes.size());

  UniformMesh &mesh = *ctx.mesh;
  UniformMesh::visitor vis, vis2;
  vis.visit_mode = UniformMesh::BreadthFirst;
  vis.init(mesh);
  vis2.init(mesh);

  int ndim = ctx.ndim;
  assert(x.size() % ndim ==0);
  ctx.nnod = x.size()/ndim;
  int nnod = ctx.nnod;

  x.reshape(2,nnod,ndim);
  int face_nel = surf_con.size(1);

  assert(ndim == 3);
  assert(face_nel == 3);
  FastMat2 X(2,face_nel,ndim), 
    edgea(1,ndim), edgeb(1,ndim), surf(1,ndim);

  int nsurf_nodes = surf_nodes.size();
  int nfaces = surf_con.size(0);
  surf_mass.a_resize(1,nfaces);
  node_mass.resize(nsurf_nodes);
  double Area = 0.;
  for (int jface=0; jface<nfaces; jface++) {
    for (int j=0; j<face_nel; j++) {
      int rnode = surf_con.e(jface,j);
      int node=surf_nodes.ref(rnode);
      assert(node<nnod);
      X.ir(1,j+1).set(&x.e(node,0));
    }
    X.rs();

    // Compute surface area for surface mass matrix
    X.ir(1,2);
    edgea.set(X);

    X.ir(1,3);
    edgeb.set(X);

    X.ir(1,1);
    edgea.rest(X);
    edgeb.rest(X);
    X.rs();

    surf.cross(edgea,edgeb);
    double area = sqrt(surf.sum_square_all())/2.0;
    surf_mass.ref(jface) = area;
    double nod_area = area/face_nel;
    for (int j=0; j<face_nel; j++) {
      int rnode = surf_con.e(jface,j);
      node_mass.ref(rnode) += nod_area;
    }
    Area += area;
    // printf("face %d, area %f\n",jface,area);
  }
  // printf("total area %f\n",Area);
}
