// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: util3.h,v 1.9 2004/04/30 15:10:46 mstorti Exp $
#ifndef PETSCFEM_UTIL3_H
#define PETSCFEM_UTIL3_H
#include <string>
#include <src/dvector.h>

#ifdef USE_SSL
#include <SSL/sockets.h>
#endif

using namespace std;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int string2int(string &s,int &n);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int string2dbl(string &s,double &b);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Converts a line in a list of tokens separated by white space. 
    #tokens# is cleared before the tokenization. 
    @param line (input) line to be tokenized
    @param tokens (output) vector f tokens
*/ 
void tokenize(const char *line,vector<string> &tokens);

// Regularized delta based on COS function, in interval [0,b]
double pf_regdelta(double x,double b);

// Regularized delta based on COS function, in interval [a,b]
double pf_regdelta(double x,double a,double b);

// Regularized Heaviside function based on COS, in interval [0,b]
double pf_regheavis(double x,double b);

// Regularized Heaviside function based on COS, in interval [a,b]
double pf_regheavis(double x,double a,double b,double y0=0.0,double y1=1.0);

double pf_regmin(double x1,double x2,double a);
double pf_regabs(double x,double a,double &dx);
double pf_regabs(double x,double a);
double pf_regabs2(double x,double a,double b);
double pf_regmin2(double y1,double y2,double a,double &mu);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define CHECK_COOKIE(keyword)							\
    { Sgetline(&buf,&Nbuf,sock);						\
    tokenize(buf,tokens);							\
    ierr = string2int(tokens[1],cookie2);					\
    PETSCFEM_ASSERT0((tokens[0]==#keyword "_OK" && !ierr && cookie==cookie2),	\
		     "Bad response from DX client sending " #keyword "\n"); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Concatenate arrays that are defined in each process
template<class T>
void concat(std::vector<T> &in,std::vector<T> &out) {
  int nproc,rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  // Compute the counts and displacements for the Allgather
  std::vector<int> counts(nproc,0), displs(nproc+1,0);
  int szT = sizeof(T);
  // We compute the sizes in bytes (this may be wrong!!)
  int nhere = in.size()*szT;
  counts[rank] = nhere;
  // We use the MPI_COMM_WORLD, I don't know if this is right
  // We gather all the local sizes so that we have in
  // each processor the sizes of all the processors
  MPI_Allgather(&nhere,1,MPI_INT,counts.data(),1,MPI_INT,MPI_COMM_WORLD);
  // Compute the displs as cumsum of the sizes
  displs[0]=0;
  for (int j=0; j<nproc; j++)
    displs[j+1] = displs[j]+counts[j];
  // Resize the output vector
  out.resize(displs[nproc]/szT);
  // Dothe Allgather
  MPI_Allgatherv(in.data(),nhere,MPI_CHAR,
                 out.data(),counts.data(),displs.data(),MPI_CHAR,
                 MPI_COMM_WORLD);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
// Concatenate arrays that are defined in each process
template<class T>
void concat(dvector<T> &in,dvector<T> &out) {
  int nproc,rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  std::vector<T> tmpin,tmpout;
  int N = in.size();
  for (int j=0; j<N; j++) 
    tmpin.push_back(in.ref(j));
  concat(tmpin,tmpout);
  if (!rank) {
    int M = tmpout.size();
    out.mono(M);
    for (int k=0; k<M; k++) out.ref(k) = tmpout[k];
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#ifdef USE_SSL
#if 0
#define SGETLINE_FACTOR 2
#define SGETLINE_INIT_SIZE 512
#define SGETLINE_MAX_SIZE INT_MAX
#endif
/** Reads a line from a socket using the Simple sockets 
    library function #Sgets# but with eventual reallocation, 
    using #malloc#. (This is similar ro the GNU #getline# function). 
    @param lineptr (input/output) the buffer where characters are read. 
    After use, you can free with #free#. 
    @param N_a (input/output) number of bytes initially allocated in #lineptr#
    @param (input) the socket where the line is read. 
    @return number of bytes read */ 
ssize_t Sgetline(char **lineptr, size_t *N_a,Socket *sock);
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class DXSplit {
  struct Subelem {
    string type;
    int subnel;
    vector<int> indices;
  };
  vector<Subelem> splitting;
public:
  /** Parses a line of the form #dx_type_1 dx_type_2 ... dx_type_n#
      where each #dx_type_j# is of the form #subel_1 subel_2 ...  dx_type#,
      #dx_type# may be #quads#, #cubes#, 
      with #k = nsubelem * subnel#. For instance
      a prism may be split in 3 tetras with a line like #tetrahedra 3 4 1 2 3 4 
      5 4 6 2 2 6 3 4#. 
      @param line (input) the line to be parsed */ 
  int parse(const char *line);

  /** Number of sub-types in the splitting. 
      @return number of subelements */ 
  int dx_types_n();

  /** Returns the description of the #j#-th type. 
      @param j (input) 0-based type index
      @param dx_type (output) the DX type
      @param subnel (input) the number of nodes for this type
      @param nodes (input) the nodes connected to this subelement. Length
      must be multiple of subnel */ 
  void dx_type(int j,string &dx_type,int &subnel,vector<int> &nodes);
  
};

#endif
