// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: dvector-epl.cpp,v 1.9 2005/05/16 03:25:05 mstorti Exp $
<:
$sfx = $type unless defined $sfx;
$dvtype="dv$sfx" unless defined $dvtype;
$macro="DV_\U$sfx" unless defined $macro;
:>//

#define <:=$macro:>
#define TYPE <:=$type:>
// #define STYPE "<:=$type:>"
#define DVTYPE "<:=$dvtype:>"

#define TAG <:=$dvtype:>_tag

#define MAKE_DVECTOR_FUN make_<:=$dvtype:>
#define FREE_DVECTOR_FUN free_<:=$dvtype:>

#define DVECTOR_CLONE_FUN <:=$dvtype:>_clone
#define DVECTOR_SCALE_FUN <:=$dvtype:>_scale_x
#ifndef DV_INT
#define DVECTOR_SCALE_FUN_S "<:=$dvtype:>-scale!"
#endif
#define DVECTOR_PUSH_FUN <:=$dvtype:>_push
#define DVECTOR_SIZE_FUN <:=$dvtype:>_size
#define DVECTOR_PRINT_FUN <:=$dvtype:>_print
#define DVECTOR_RESIZE_FUN <:=$dvtype:>_resize
#define DVECTOR_SET_W1_FUN <:=$dvtype:>_set_w1
#define DVECTOR_SET_W2_FUN <:=$dvtype:>_set_w2
#define DVECTOR_READ_FUN <:=$dvtype:>_read_x
#define DVECTOR_CAT_FUN <:=$dvtype:>_cat_x
#define DVECTOR_DUMP_FUN <:=$dvtype:>_dump1
#define DVECTOR_RESHAPE_FUN <:=$dvtype:>_reshape
#define DVECTOR_SHAPE_FUN <:=$dvtype:>_shape

#define INIT_DVECTOR_FUN <:=$dvtype:>_init

#include "./dvector.cpp"
