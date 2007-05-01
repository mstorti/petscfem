# -*- mode: C++ -*-
# $Id: fmat2.pl,v 1.2 2003/07/02 03:36:13 mstorti Exp $
$cache_op=<<'//EOF';
FastMatCache *cache = ctx->step();
//EOF

$CACHE_OPERATIONS = $cache_op;
