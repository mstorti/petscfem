// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: fifo.h,v 1.1 2002/12/30 03:06:31 mstorti Exp $
#ifndef FIFO_H
#define FIFO_H

#define _GNU_SOURCE

#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>

void read_doubles(const char *line,const char *name,
		  vector<double> &v);

void read_doubles(FILE *fid,const char *name, vector<double> &v,
		  int n=-1);

double read_doubles(FILE *fid,const char *name);

void read_doubles2(FILE *fid,const char *name, vector<double> &v,
		  int n=-1);

double read_doubles2(FILE *fid,const char *name);

const char* read_string(FILE *fid);

#endif
