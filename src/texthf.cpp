#include "texthf.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTableFilter::get_entry"
void TextHashTableFilter::get_entry(const char *kk,const char *&val) {
  string key;
  vector<string>::iterator k;
  for (k=prefix.end(); k!=prefix.begin(); k--)
    key += *k + string(".");
  key += string(kk);
  get_entry(key.c_str(),val);
}
