// -*- mode: c++ -*-
//__INSERT_LICENSE__
//$Id: inviscid.h,v 1.4 2003/01/04 00:49:56 mstorti Exp $
#ifndef ROSI_H
#define ROSI_H

#include <src/ampli.h>
#include <cstdio>
#include <sys/stat.h>
#include <cassert>
#include <cstdlib>
#include <vector>

#include <src/vecmacros.h>
#include <src/fstack.h>
#include <src/texthash.h>
#include <src/fem.h>
#include <src/texthf.h>
#include <src/sttfilter.h>

#include <applications/ns/nsi_tet.h>
#include <applications/ns/nssup.h>
#include <applications/ns/dlhook.h>

//#include "./fifo.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class coupling_visc_hook {
private:
  FILE *visc2inv, *inv2visc;
public:
  void init(Mesh &mesh,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close();
};

vector<double> xnod, u;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class coupling_inv_hook {
private:
  int nnod;
  FILE *visc2inv, *inv2visc;
public:
  void init(Mesh &mesh,Dofmap &dofmap,
	    TextHashTableFilter *options,const char *name);
  void time_step_pre(double time,int step);
  void time_step_post(double time,int step,
		      const vector<double> &gather_values);
  void close() {}
};

#endif
