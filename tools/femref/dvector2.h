// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: dvector2.h,v 1.1 2005/01/17 02:47:47 mstorti Exp $
#ifndef SCHEME_FEM_DVECTOR2_H
#define SCHEME_FEM_DVECTOR2_H

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <libguile.h>

#include <src/dvector.h>

void 
scmlist2vec(SCM s_list,vector<int> &v);

#endif
