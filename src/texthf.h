#ifndef TEXTHF_H
#define TEXTHF_H

#include <map>
#include <vector>
#include <string>
#include <cstring>
#include "getprop.h"
#include "texthash.h"

class TextHashTableFilter {
private:
  TextHashTable * thash;
  vector<string> prefix;
  mutable string key,result;
  void set_key(const char *name,string &s) const;
public:
  TextHashTableFilter(TextHashTable *thash_p) : thash(thash_p) {};
  TextHashTableFilter &push(const char *s) {prefix.push_back(string(s)); return *this;};
  TextHashTableFilter &pop() {prefix.pop_back(); return *this;};
  TextHashTableFilter &clear() {prefix.clear(); return *this;}
  void get_entry(const char *,const char *&) const;
  int get_int(const char *name,
	      int &retval,int defval=0,int n=1) const {
    set_key(name,key);
    return ::get_int(thash,key.c_str(),&retval,defval,n);
  };
  int get_double(const char *name,
		 double &retval,int defval=0,int n=1) const {
    set_key(name,key);
    return ::get_double(thash,key.c_str(),&retval,defval,n);
  };
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets a vector of double properties whose length is unknown. 
      @author M. Storti
      @param name (input) the name of the property
      @param retval (input) the returned vector
      @param defval (input) Controls the action if no entry is found
      in the hash table. If #defval==0#, let retval unchanged, so that
      you have to set it to its default value. Give an error if
      #defval!=0# 
  */ 
#if 0
  int get_vec_double(const char *name,
		     vector<double> &retval,int defval=0) const;
#endif
  int get_string(const char *name,
		 string &ret,int defval=0,int n=1) const {
    set_key(name,key);
    return ::get_string(thash,key.c_str(),ret,defval,n);
  };
};

#endif
