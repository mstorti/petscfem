/*__INSERT_LICENSE__*/
//$Id: thash2.cpp,v 1.1 2002/11/03 09:18:20 mstorti Exp $

#include <cstdio>
#include <src/fstack.h>
#include <src/texthash.h>
#include <src/utils.h>

int main () {
  char *line;
  TextHashTable *t;
  FileStack *f = new FileStack("thash2.txt");
  while(f->get_line(line)) {
    if (strcmp(line,"table t1")) break;
  }
  t = new TextHashTable;
  t->read(f);
  t->print("t1:");
  f->close();
  delete f;
}
