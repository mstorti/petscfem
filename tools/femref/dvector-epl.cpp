// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: dvector-epl.cpp,v 1.5 2005/01/17 02:47:47 mstorti Exp $
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

#define DVECTOR_PUSH_FUN <:=$dvtype:>_push
#define DVECTOR_SIZE_FUN <:=$dvtype:>_size
#define DVECTOR_PRINT_FUN <:=$dvtype:>_print
#define DVECTOR_RESIZE_FUN <:=$dvtype:>_resize
#define DVECTOR_SET_FUN <:=$dvtype:>_set
#define DVECTOR_READ_FUN <:=$dvtype:>_read_x
#define DVECTOR_CAT_FUN <:=$dvtype:>_cat_x
#define DVECTOR_DUMP_FUN <:=$dvtype:>_dump
#define DVECTOR_RESHAPE_FUN <:=$dvtype:>_reshape

#define INIT_DVECTOR_FUN <:=$dvtype:>_init

#include "./dvector.cpp"
