/*__INSERT_LICENSE__*/
// $Id: tryme.cpp,v 1.1 2001/11/28 02:17:31 mstorti Exp $

#include <cassert>
#include <cstdio>
#include <pthread.h>

void *thread_main(void *p) {
  int &result = *(int *)p;
  result = 0;
  for (int k=0; k<10000; k++) result += k;
  return NULL;
}

int main() {
  pthread_t thread;
  int result;
  void *retval;
  int ierr = pthread_create(&thread,NULL,&thread_main,&result);
  assert(!ierr);
  ierr = pthread_join(thread,&retval);
  printf("result %d\n",result);
}
