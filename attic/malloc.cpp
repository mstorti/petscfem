// This was wased in NS for detecting a bug in
// the memory management

//#define DEBUG_MALLOC_USE
#ifdef DEBUG_MALLOC_USE

#include <map>
map<void *,unsigned int> mem_map;

FILE* malloc_log;
unsigned int total=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/* Global variables used to hold underlaying hook values.  */
static void *(*old_malloc_hook) (size_t,const void *);
static void (*old_free_hook) (void*,const void*);
     
/* Prototypes for our hooks.  */
static void *my_malloc_hook (size_t,const void*);
static void my_free_hook(void*,const void*);
     
static void *
my_malloc_hook (size_t size,const void *caller)
{
  void *result;
  static unsigned int last=0,chunk_size=1000000;
  /* Restore all old hooks */
  __malloc_hook = old_malloc_hook;
  __free_hook = old_free_hook;
  /* Call recursively */
  result = malloc (size);
  /* Save underlaying hooks */
  old_malloc_hook = __malloc_hook;
  old_free_hook = __free_hook;
  /* `printf' might call `malloc', so protect it too. */
  //  fprintf (malloc_log,"malloc (%u) returns %p\n", (unsigned int)
  // size, result);
  mem_map[result] = (unsigned int) size;
  total += size;
  if (total >(last+1)*chunk_size) {
    printf("allocated %d\n",total);
    last++;
  }
  /* Restore our own hooks */
  __malloc_hook = my_malloc_hook;
  __free_hook = my_free_hook;
  return result;
}
  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:    
static void my_free_hook (void *ptr,const void *caller)
{
  /* Restore all old hooks */
  __malloc_hook = old_malloc_hook;
  __free_hook = old_free_hook;
  /* Call recursively */
  free (ptr);
  /* Save underlaying hooks */
  old_malloc_hook = __malloc_hook;
  old_free_hook = __free_hook;
  /* `printf' might call `free', so protect it too. */
  //  fprintf (malloc_log,"freed pointer %p\n", ptr);
  if (mem_map.find(ptr) == mem_map.end()) {
    printf("Not found this pointer!! %p\n",ptr);
  } else {
    unsigned int s = mem_map[ptr];
    mem_map.erase(ptr);
    total -= s;
  }
  /* Restore our own hooks */
  __malloc_hook = my_malloc_hook;
  __free_hook = my_free_hook;
}
#endif
