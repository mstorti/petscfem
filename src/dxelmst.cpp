//__INSERT_LICENSE__
//$Id: dxelmst.cpp,v 1.1 2003/02/11 21:37:00 mstorti Exp $

#include <vector>

#include <src/fem.h>
#include <src/elemset.h>
#include <src/util3.h>

#ifdef USE_DX
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::dx(Socket *sock,Nodedata *nd,double *field_state) {
  int ierr, cookie, cookie2;

#define BUFSIZE 512
  static char *buf = (char *)malloc(BUFSIZE);
  static size_t Nbuf = BUFSIZE;
  vector<string> tokens;

  //o Flags whether the element should return a connection array
  //   for this elemset. 
  TGETOPTDEF(thash,int,dx,0);
  if (!dx) return;

  string type;
  int subnel,nsubelem;
  vector<int> node_indices;
  assert(dx_types_n()==1);

  dx_type(0,type,subnel,node_indices);
  assert(subnel>0);
  assert(node_indices.size() % subnel == 0);
  nsubelem = node_indices.size()/subnel;

  if (!MY_RANK) {
    cookie = rand();
    Sprintf(sock,"elemset %s %s %d %d %d\n",name(),type.c_str(),
	    subnel,nelem*nsubelem,cookie);
    printf("Sending elemset %s %s %d %d %d\n",name(),type.c_str(),
	   subnel,nelem*nsubelem,cookie);
    for (int j=0; j<nelem; j++) {
      int *row = icone+j*nel;
      for (int jj=0; jj<nsubelem; jj++) {
	for (int n=0; n<subnel; n++) {
	  int k = node_indices[jj*subnel+n];
	  // Convert to 0 based (DX) node numbering
	  int node = *(row+k)-1;
	  Swrite(sock,&node,sizeof(int));
	}
      }
    }
    CHECK_COOKIE(elemset);
    cookie = rand();
    Sprintf(sock,"field %s_field nodes %s state %d\n",name(),name(),cookie);
    CHECK_COOKIE(field);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int Elemset::dx_types_n() {
  const char *line;
  //o _T: <string:type> <integer array>
  //  _N: dx_indices 
  //  _D: (empty array)
  //  _DOC: list of indices of the nodes to be passed to DX 
  //        as connectivity table. Note that the order the nodes
  //        are entered in DX is not the same that in PETSc-FEM.
  //        For instance, in  PETSc-FEM the nodes of a quad are entered
  //        counter-clockwise, and it they are (in that order) 1, 2, 3, 4, then
  //        for DX they are entered in the order: 1, 2, 4, 3. Also DX wants 0-based
  //        (C-sytle) node numbers, whereas PETSc-FEM uses 1-based (Fortran style).
  //        However they must be entered 1-based for this option. For instance,
  //        for quads one should enter {\tt dx_indices 1 2 4 3} and for cubes
  //        {\tt dx_indices 1 2 4 3 5 6 8 7}. 
  //  _END
  thash->get_entry("dx_indices",line);
  // Reads list of indices from `dx_indices' option line
  // othewise, use standard numeration from geometry
  if (!line) {
    assert(0);
    // for (int k=0; k<nel; k++) node_indices.push_back(k);
  } else {
    splitting.parse(line);
    return splitting.dx_types_n();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::dx_type(int j,string &dx_type,int &subnel,vector<int> &nodes) {
  return splitting.dx_type(j,dx_type,subnel,nodes);
}

#endif
