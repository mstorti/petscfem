//__INSERT_LICENSE__
// $Id: getsurf.cpp,v 1.16 2005/01/14 23:26:26 mstorti Exp $

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
#include <libguile.h>

using namespace std;

#define VERBOSE 0

#include "./femref.h"
#include "./gtemplates.h"

typedef SCM(*scm_fun)();

struct FaceIterator {
  int elem;
  int face;
};

typedef multimap<int,FaceIterator> face_table_t;
typedef pair<int,FaceIterator> ft_pair_t;

extern "C"
void getsurf(SCM s_iconef,SCM s_xnodf, SCM s_statef,
	     SCM s_sconf, SCM s_graduf, SCM s_base) {

  const char *iconef = SCM_STRING_CHARS(s_iconef);
  const char *xnodf  = SCM_STRING_CHARS(s_xnodf);
  const char *statef = SCM_STRING_CHARS(s_statef);
  const char *sconf = SCM_STRING_CHARS(s_sconf);
  const char *graduf = SCM_STRING_CHARS(s_graduf);
  int base = SCM_INUM(s_base);
  int nread; 
  char c;

  assert(xnodf);
  assert(iconef);
  int ndof = 4, ndim = 3;
  time_t start, end;
  dvector<double> u,x;
  const GeomObject::Template *mesh_tmpl = &OrientedTetraTemplate;
  UniformMesh mesh(*mesh_tmpl,3);
  int mesh_nel = mesh.tmplt()->size_m;
  mesh.read(xnodf,iconef,base);
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
    X(2,mesh_nel,ndim+1), invX,
    edgea(1,ndim), edgeb(1,ndim), surf(1,ndim);
  X.ir(2,ndim+1).set(1.0).rs();
  FILE *fid = fopen(sconf,"w");
  FILE *fidgu = fopen(graduf,"w");
  map<int,int> surf_nodes_map;
  int nfaces = face_table.size(),
    nsurf_nodes = 0;
  dvector<double> grad_Ue;
  grad_Ue.a_resize(2,nfaces,ndim*ndof);
  grad_Ue.defrag();
  dvector<int> surf_nodes;
  dvector<int> surf_con;
  dvector<double> surf_mass, node_mass;
  surf_con.a_resize(2,nfaces,face_nel);
  surf_mass.a_resize(1,nfaces);
  int surf_elem=0;
  double Area = 0.;
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
    grad_U.export_vals(&grad_Ue.e(surf_elem,0));
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

    // Compute surface area for surface mass matrix
    X.is(2,1,ndim).ir(1,2);
    edgea.set(X);

    X.ir(1,3);
    edgeb.set(X);

    X.ir(1,1);
    edgea.rest(X);
    edgeb.rest(X);
    X.rs();

    surf.cross(edgea,edgeb);
    double area = sqrt(surf.sum_square_all())/2.0;
    surf_mass.ref(surf_elem) = area;
    Area += area;
    surf_elem++;
  }
  surf_nodes.defrag();
  fclose(fid);
  fclose(fidgu);
  printf("total area %f\n",Area);
  node_mass.resize(nsurf_nodes);
  dvector<double> grad_Un;
  grad_Un.a_resize(2,nsurf_nodes,ndim*ndof);
  grad_Un.defrag();

  // Save surface connectivity (reduced representation)
  fid = fopen("cube.surf-con-red.tmp","w");
  for (int jface=0; jface<nfaces; jface++) {
    for (int j=0; j<face_nel; j++) {
      int node = surf_con.e(jface,j);
      fprintf(fid,"%d ",surf_con.e(jface,j));
    }
    fprintf(fid,"\n");
  }
  fclose(fid);

  // Save surface node coordinates
  fid = fopen("cube.surf-nod.tmp","w");
  for (int jnode=0; jnode<nsurf_nodes; jnode++) {
    int node = surf_nodes.ref(jnode);
    for (int j=0; j<ndim; j++) 
      fprintf(fid,"%g ",x.e(node,j));
    fprintf(fid,"\n");
  }
  fclose(fid);

  // Loop over smoothing steps
  int niter=2, jiter=0;
  while (true) {
    grad_Un.set(0.);
    for (int jface=0; jface<nfaces; jface++) {
      double nod_area = surf_mass.ref(jface)/double(face_nel);
      for (int j=0; j<face_nel; j++) {
	int node = surf_con.e(jface,j);
	if (jiter==0) node_mass.ref(node) += nod_area;
	double 
	  *to = &grad_Un.e(node,0),
	  *from = &grad_Ue.e(jface,0);
	for (int k=0; k<ndim*ndof; k++) 
	  to[k] += from[k]*nod_area;
      }
    }
    for (int node=0; node<nsurf_nodes; node++) {
      double nod_area = node_mass.ref(node);
      for (int k=0; k<ndim*ndof; k++) 
	grad_Un.e(node,k) /= nod_area;
    }
    jiter++;
    if (jiter==niter) break;

    grad_Ue.set(0.);
    for (int jface=0; jface<nfaces; jface++) {
      double *to = &grad_Ue.e(jface,0);
      for (int j=0; j<face_nel; j++) {
	int node = surf_con.e(jface,j);
	double *from = &grad_Un.e(node,0);
	for (int k=0; k<ndim*ndof; k++) 
	  to[k] += from[k];
      }
      for (int k=0; k<ndim*ndof; k++) 
	to[k] /= double(face_nel);
    }
  }

  fid = fopen("cube.grad-un.tmp","w");
  for (int node=0; node<nsurf_nodes; node++) {
    // printf("%d ",node);
    for (int k=0; k<ndim*ndof; k++) 
      fprintf(fid,"%g ",grad_Un.e(node,k));
    fprintf(fid,"\n");
  }
  fclose(fid);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
extern "C" void
init_femref(void) {
  scm_c_define_gsubr("getsurf",6,0,0,scm_fun(getsurf));
}
