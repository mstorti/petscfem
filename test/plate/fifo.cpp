/*__INSERT_LICENSE__*/
// $Id: fifo.cpp,v 1.1 2002/12/30 03:06:31 mstorti Exp $
#include "./fifo.h"
#include <cstdlib>

void read_doubles(const char *line,const char *name,
		  vector<double> &v) {
  char *tok;
  double val;
  int nread;

  v.clear();
  char *line_copy = (char *)malloc((strlen(line)+1)*sizeof(char));
  strcpy(line_copy,line);
  tok = strtok(line_copy," ");
  if (!tok || strcmp(tok,name)) {
    printf("read_doubles: bad line or line doesn't match key\n"
	   "line: \"%s\"\n"
	   "key: \"%s\"\n",line,name);
  }
  tok = strtok(NULL," ");
  while (tok) {
    nread = sscanf(tok,"%lf",&val);
    assert(nread<=1);
    if (nread<0) break;
    v.push_back(val);
    tok = strtok(NULL," ");
  }
}

void read_doubles(FILE *fid,const char *name, vector<double> &v,
		  int n=-1) {
  static char *line=NULL;
  static size_t N=0;
  // if (!line) 
  // line = (char *)malloc(N*sizeof(char));
  int nread = getline(&line, &N, fid);
  assert(nread!=-1);
  read_doubles(line,name,v);
  if (n>=0) assert(v.size()==n);
}

double read_doubles(FILE *fid,const char *name) {
  vector<double> v;
  read_doubles(fid,name,v,1);
  return v[0];
}

double read_doubles2(FILE *fid,const char *name) {
  vector<double> v;
  read_doubles2(fid,name,v,1);
  return 0.;
}

void read_doubles2(FILE *fid,const char *name, vector<double> &v,
		  int n=-1) {
  static char *line=NULL;
  static size_t N=0;
  // if (!line) 
  // line = (char *)malloc(N*sizeof(char));
  int nread = getline(&line, &N, fid);
  assert(nread!=-1);
  printf("wants to read key \"%s\", reads \"%s\"\n",name,line);
  if (n>0) v.resize(n,0.);
}

const char* read_string(FILE *fid) {
  static char *line=NULL;
  static size_t N=0;
  int nread = getline(&line, &N, fid);
  assert(nread!=-1);
  // Chomp end-of-line
  for (int j=0; j<strlen(line); j++) {
    if (line[j]=='\n') {
      line[j] = '\0';
      break;
    }
  }
  return line;
}

