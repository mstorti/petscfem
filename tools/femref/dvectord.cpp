// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: dvectord.cpp,v 1.10 2005/01/17 15:43:26 mstorti Exp $

#define DV_DBL
#define TYPE double
// #define STYPE "double"
#define DVTYPE "dvdbl"

#define TAG dvdbl_tag

#define MAKE_DVECTOR_FUN make_dvdbl
#define FREE_DVECTOR_FUN free_dvdbl

#define DVECTOR_PUSH_FUN dvdbl_push
#define DVECTOR_SIZE_FUN dvdbl_size
#define DVECTOR_PRINT_FUN dvdbl_print
#define DVECTOR_RESIZE_FUN dvdbl_resize
#define DVECTOR_SET_W1_FUN dvdbl_set_w1
#define DVECTOR_SET_W2_FUN dvdbl_set_w2
#define DVECTOR_READ_FUN dvdbl_read_x
#define DVECTOR_CAT_FUN dvdbl_cat_x
#define DVECTOR_DUMP_FUN dvdbl_dump
#define DVECTOR_RESHAPE_FUN dvdbl_reshape
#define DVECTOR_SHAPE_FUN dvdbl_shape

#define INIT_DVECTOR_FUN dvdbl_init

#include "./dvector.cpp"
