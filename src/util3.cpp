//__INSERT_LICENSE__
// $Id: util3.cpp,v 1.3 2003/02/09 14:50:57 mstorti Exp $
#include <cstring>
#include <cstdio>
#include <string>
#include <vector>
#include <mpi.h>

#include <HDR/sockets.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int string2int(string s,int &n) {
  int nread = sscanf(s.c_str(),"%d",&n);
  return (nread!=1);
} 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Converts a line in a list of tokens separated by white space. 
    #tokens# is cleared before the tokenization. 
    @param line (input) line to be tokenized
    @param tokens (output) vector f tokens
*/ 
void tokenize(const char *line,vector<string> &tokens) {
  // Make a local copy (input is read only)
  char *copy = new char[strlen(line)+1];
  strcpy(copy,line);
  // White space pattern
  char spc[] = "[ \t\n]";
  // Clear tokens arg
  tokens.clear();
  // Tokenize using `strtok'
  int j=0;
  while(1) {
    char *token = strtok((j ? NULL : copy),spc);
    if (!token) break;
    tokens.push_back(token);
    j++;
  }
  // clear local copy
  delete[] copy;
}

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
ssize_t Sgetline(char **lineptr, size_t *N_a,Socket *sock) {
  unsigned int &N = *N_a;	// use reference for better readbility
  char * new_line_ptr = NULL, *q, *q0, *qe;
  // At any time the buffer is #N# bytes long and we have
  // read already #read_so_far# bytes.
  int read_so_far=0;
  // Main loop. We read lines with gets until a "\n" is found. If the line
  // has not a "\n" then it should end in "\0\0". 
  while (1) {
    if (N>0) {
      // We read onto lineptr[q0,qe)
      q0 = *lineptr+read_so_far;
      qe = *lineptr+N;
      // Set all to nulls (if a `\n' is left, then we could detect
      // a false line. 
      for (q = q0; q< qe; q++) *q = '\0';
      // Get next part of the line
      Sgets(q0,N-read_so_far,sock);
      // pprint(*lineptr,N);
      // If a newline is found, then we have read the line
      for (q = q0; q<qe; q++) if (*q == '\n') break;
      if (q<qe) break;
      // If not end, then we have read all, except the two nulls at the end
      read_so_far = N-2;
    }
    // Realocate
    N = (N ? 2*N : SGETLINE_INIT_SIZE);
    if (N > SGETLINE_MAX_SIZE) return 0;
    new_line_ptr = (char *) malloc(N);
    assert(new_line_ptr);
    // If already have a buffer, copy to new allocated. 
    if (*lineptr) {
      memcpy(new_line_ptr,*lineptr,read_so_far);
      free(*lineptr);
    }
    // update pointer
    *lineptr = new_line_ptr;
    new_line_ptr = NULL;
  }
  // return number of bytes read
  return strlen(*lineptr)+1;
}

