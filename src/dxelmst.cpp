//__INSERT_LICENSE__
//$Id: dxelmst.cpp,v 1.7 2003/09/10 23:18:43 mstorti Exp $

#ifdef USE_DX
#include <vector>

#include <src/fem.h>
#include <src/elemset.h>
#include <src/util3.h>
#include <src/sockbuff.h>

extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Elemset::dx_send_connectivities(Socket *sock,
				     int nsubelem, int subnel,
				     vector<int> &node_indices) {
  SocketBuffer<int> sbuff(sock);
  for (int j=0; j<nelem; j++) {
    int *row = icone+j*nel;
    for (int jj=0; jj<nsubelem; jj++) {
      for (int n=0; n<subnel; n++) {
	int k = node_indices[jj*subnel+n];
	// Convert to 0 based (DX) node numbering
	sbuff.put(*(row+k)-1);
	    }
    }
  } 
  sbuff.flush();
}

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

  //o Uses DX cache if possible in order to avoid sending the connectivities
  //  each time step or frame. Use only if the connectivities are not changing
  //  in your problem. 
  TGETOPTDEF(thash,int,dx_cache_connectivities,0);

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
    if (dx_cache_connectivities) {
      Sprintf(sock,"elemset %s %s %d %d %d use_cache\n",name(),type.c_str(),
	      subnel,nelem*nsubelem,cookie);
      
      Sgetline(&buf,&Nbuf,sock);
      tokenize(buf,tokens);
      
      if (tokens[0]=="send_elemset") {
	printf("Sending elemset...\n");
	dx_send_connectivities(sock,nsubelem,subnel,node_indices);
      } else if (tokens[0]=="do_not_send_elemset") {
	printf("Does not send elemset.\n");
      } else PETSCFEM_ERROR("Error in DXHOOK protocol. DX sent \"%s\"\n",
			    tokens[0].c_str());
    } else {
      Sprintf(sock,"elemset %s %s %d %d %d\n",name(),type.c_str(),
	      subnel,nelem*nsubelem,cookie);
      dx_send_connectivities(sock,nsubelem,subnel,node_indices);
    }
    CHECK_COOKIE(elemset);
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
