// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: hashf.h,v 1.1 2004/01/21 02:47:37 mstorti Exp $
#ifndef PETSCFEM_HASHF_H
#define PETSCFEM_HASHF_H

typedef  unsigned long  int  ub4;   /* unsigned 4-byte quantities */
typedef  unsigned       char ub1;   /* unsigned 1-byte quantities */

#define hashsize(n) ((ub4)1<<(n))
#define hashmask(n) (hashsize(n)-1)
ub4 hash_fun(ub1 *k,ub4 length,ub4 initval);

#endif
