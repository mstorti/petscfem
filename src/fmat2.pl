# -*- mode: C++ -*-
# $Id: fmat2.pl,v 1.1 2002/11/28 15:13:25 mstorti Exp $
$cache_op=<<'//EOF';
FastMatCache *cache;

if (was_cached) {
  cache = cache_list_begin[position_in_cache++];
#ifdef FM2_CACHE_DBG
  printf ("reusing cache: ");
#endif
} else if (!use_cache) {
  cache = new FastMatCache;
} else {
  cache = new FastMatCache;
  cache_list->push_back(cache);
  cache_list_begin = cache_list->begin();
  cache_list->list_size =
    cache_list_size = cache_list->size();
  position_in_cache++;
#ifdef FM2_CACHE_DBG
  printf ("defining cache: ");
#endif
}
#ifdef FM2_CACHE_DBG
printf(" cache_list %p, cache %p, position_in_cache %d\n",
       cache_list,cache,position_in_cache-1);
#endif
//EOF

$CACHE_OPERATIONS = $cache_op;
