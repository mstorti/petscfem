// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: util3.h,v 1.3 2003/02/09 14:50:57 mstorti Exp $
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

#endif
