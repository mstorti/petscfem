//__INSERT_LICENSE__
//$Id: dxelmst.cpp,v 1.2 2003/02/12 00:36:27 mstorti Exp $

#ifdef USE_DX
#include <vector>

#include <src/fem.h>
#include <src/elemset.h>
#include <src/util3.h>

extern int MY_RANK,SIZE;

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
  int ierr;
  //o _T: mixed string/integer array
  //  _N: dx_indices 
  //  _D: (none)
  //  _DOC: 
  //i_tex ../doc/nsdoc.tex dx_line
  //  _END
  thash->get_entry("dx_indices",line);
  // Reads list of indices from `dx_indices' option line
  // othewise, use standard numeration from geometry
  if (!line) {
    // Uses splitting from the GPdata object
    const char *geom;
    thash->get_entry("geometry",geom);
    assert(geom);
    // Number of Gauss points.
    TGETOPTDEF(thash,int,npg,0); //nd
    // Dimension of the problem
    TGETOPTDEF(thash,int,ndim,0); //nd
    // Dimension of the element
    TGETOPTDEF(thash,int,ndimel,0); //nd
    if (ndimel==0) ndimel=ndim;
    GPdata gpdata(geom,ndimel,nel,npg);
    splitting = gpdata.splitting;
  } else {
    splitting.parse(line);
  }
  return splitting.dx_types_n();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::dx_type(int j,string &dx_type,int &subnel,vector<int> &nodes) {
  return splitting.dx_type(j,dx_type,subnel,nodes);
}

#endif
