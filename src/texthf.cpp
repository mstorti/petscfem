//__INSERT_LICENSE__
//$Id: texthf.cpp,v 1.4 2003/01/08 15:54:25 mstorti Exp $
#include "texthf.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTableFilter::set_key"
void TextHashTableFilter::set_key(const char *name,string &kk) const {
  vector<string>::const_iterator k;
  kk.erase();
  for (k=prefix.begin(); k!=prefix.end(); k++) {
    kk +=  *k; 
    kk += string(".");
  }
  kk += string(name);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTableFilter::get_entry"
void TextHashTableFilter::get_entry(const char *kk,const char *&val)
  const {
  set_key(kk,key);
  thash->get_entry(key.c_str(),val);
}
