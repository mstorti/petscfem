// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: fifo.h,v 1.2 2003/07/05 03:28:35 mstorti Exp $
#ifndef FIFO_H
#define FIFO_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>

using namespace std;

void read_doubles(const char *line,const char *name,
		  vector<double> &v);

void read_doubles(FILE *fid,const char *name, vector<double> &v,
		  int n);

double read_doubles(FILE *fid,const char *name);

void read_doubles2(FILE *fid,const char *name, vector<double> &v,
		  int n=-1);

double read_doubles2(FILE *fid,const char *name);

const char* read_string(FILE *fid);

#endif
