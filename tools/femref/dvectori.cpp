// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: dvectori.cpp,v 1.1 2005/01/15 11:55:48 mstorti Exp $

#define DV_INT
#define TYPE int
// #define STYPE "int"
#define DVTYPE "dvint"

#define TAG dvint_tag

#define MAKE_DVECTOR_FUN make_dvint
#define FREE_DVECTOR_FUN free_dvint

#define DVECTOR_PUSH_FUN dvint_push
#define DVECTOR_SIZE_FUN dvint_size
#define DVECTOR_PRINT_FUN dvint_print
#define DVECTOR_RESIZE_FUN dvint_resize
#define DVECTOR_SET_FUN dvint_set
#define DVECTOR_REF_FUN dvint_ref

#define INIT_DVECTOR_FUN dvint_init

#include "./dvector.cpp"
