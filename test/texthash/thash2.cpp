/*__INSERT_LICENSE__*/
//$Id: thash2.cpp,v 1.2 2002/11/03 09:36:53 mstorti Exp $

#include <cstdio>
#include <src/fstack.h>
#include <src/texthash.h>
#include <src/utils.h>

int main () {
  char *line;
  TextHashTable *t;

  // Reads tables from file
  FileStack *f = new FileStack("thash2.txt");
  while(!f->get_line(line)) {
    if (strcmp(line,"table t1")) break;
  }
  t = new TextHashTable;
  t->read(f);
  t->print("t1:");
  t->register_name("t1");

  while(!f->get_line(line)) {
    if (strcmp(line,"table t2")) break;
  }
  t = new TextHashTable;
  t->read(f);
  t->print("t2:");
  t->register_name("t2");

  f->close();
  delete f;

  // Look for table t2 and print
  const TextHashTable *tt = TextHashTable::find("t2");
  tt->print("should be t2:");

  // Look for table t1 and print
  tt = TextHashTable::find("t1");
  tt->print("should be t1:");
}
