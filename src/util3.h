// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: util3.h,v 1.7 2003/07/02 23:22:19 mstorti Exp $
#ifndef PETSCFEM_UTIL3_H
#define PETSCFEM_UTIL3_H
#include <string>
#ifdef USE_SSL
#include <HDR/sockets.h>
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define CHECK_COOKIE(keyword)							\
    { Sgetline(&buf,&Nbuf,sock);						\
    tokenize(buf,tokens);							\
    ierr = string2int(tokens[1],cookie2);					\
    PETSCFEM_ASSERT0((tokens[0]==#keyword "_OK" && !ierr && cookie==cookie2),	\
		     "Bad response from DX client sending " #keyword "\n"); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#ifdef USE_SSL
#define SGETLINE_FACTOR 2
#define SGETLINE_INIT_SIZE 512
#define SGETLINE_MAX_SIZE INT_MAX
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
      a prism may be slit in 3 tetras with a line like #tetrahedra 3 4 1 2 3 4 
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
