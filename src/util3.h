// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: util3.h,v 1.4 2003/02/11 11:34:00 mstorti Exp $
#ifndef UTIL3_H
#define UTIL3_H

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int string2int(string s,int &n);

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
      where each #dx_type_j# is of the form #dx_type nsubelem subnel n1 n2 ... nk#,
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
