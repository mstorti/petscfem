#define _GNU_SOURCE
#include <iostream>

#include <cassert>
#include <cmath>
#include <string>

int main() {
  string s;
  s << 23 << "  " << &s;
  cout << 
}

#if 0
const char *rand_mes() {
  string s;
  static char *mes = NULL;
  // if (mes==NULL) mes = (char *)malloc(sizeof(char)*20);
  int len = irand(1,100);
  s = string("");
  for (int j=0; j<len; j++) 
    s = s + char('a'+j % 26);
  free(mes);
  int ret = asprintf(&mes,"string \"%s\"",s.c_str());
  assert(ret>=0);
  return mes;
}

int main() {
//    for (int j=0; j<20; j++)
//      printf("%s\n",rand_mes());
  const char *s;
  int j=0;
  while(1) {
    s = rand_mes();
    if (++j % 1000 ==0 ) printf("%s\n",s);
  }
}
#endif
