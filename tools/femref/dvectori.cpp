// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: dvectori.cpp,v 1.16 2005/03/03 02:20:56 mstorti Exp $

#define DV_INT
#define TYPE int
// #define STYPE "int"
#define DVTYPE "dvint"

#define TAG dvint_tag

#define MAKE_DVECTOR_FUN make_dvint
#define FREE_DVECTOR_FUN free_dvint

#define DVECTOR_CLONE_FUN dvint_clone
#define DVECTOR_SCALE_FUN dvint_scale_x
#ifndef DV_INT
#define DVECTOR_SCALE_FUN_S "dvint-scale!"
#endif
#define DVECTOR_PUSH_FUN dvint_push
#define DVECTOR_SIZE_FUN dvint_size
#define DVECTOR_PRINT_FUN dvint_print
#define DVECTOR_RESIZE_FUN dvint_resize
#define DVECTOR_SET_W1_FUN dvint_set_w1
#define DVECTOR_SET_W2_FUN dvint_set_w2
#define DVECTOR_READ_FUN dvint_read_x
#define DVECTOR_CAT_FUN dvint_cat_x
#define DVECTOR_DUMP_FUN dvint_dump
#define DVECTOR_RESHAPE_FUN dvint_reshape
#define DVECTOR_SHAPE_FUN dvint_shape

#define INIT_DVECTOR_FUN dvint_init

#include "./dvector.cpp"
