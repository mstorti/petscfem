// -*- mode: C++ -*-
//__INSERT_LICENSE__
// $Id: dvector-epl.cpp,v 1.1 2005/01/15 11:55:48 mstorti Exp $
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
#define DVECTOR_REF_FUN <:=$dvtype:>_ref

#define INIT_DVECTOR_FUN <:=$dvtype:>_init

#include "./dvector.cpp"
