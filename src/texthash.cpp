/*__INSERT_LICENSE__*/
//$Id: texthash.cpp,v 1.7 2001/05/30 03:58:50 mstorti Exp $
 
//  #include <stdio.h>
//  #include <string.h>
//  #include "util2.h"
#include "texthash.h"

// The static member(s)
THashTable TextHashTable::thash_table;
TextHashTable *TextHashTable::global_options=NULL;

int TextHashTable::print_statistics=0;


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTableVal::TextHashTableVal(char *)"
TextHashTableVal::TextHashTableVal(char *s_=NULL) : called_times(0) {
  if (!s_) {
    s=s_;
  } else {
    s=local_copy(s_);
  }
};


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void print_hash_entry(void *p, void *q, void *u)"
void print_hash_entry(void *p, void *q, void *count) {
  char *pp; TextHashTableVal *qq;
  pp= (char *)p;
  qq= (TextHashTableVal *)q;
  printf("%20s -> %-10s",pp,qq->s);
  if (TextHashTable::print_statistics) {
    int n=qq->called_times;
    printf(" [%d%s]",n,(n==0 ? " - NOT ACCESSED!!" : 
			  n<INT_MAX ? "" : "+"));
  }
  printf("\n");
  (*(int *)count)++;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void delete_hash_entry(void *p, void *q, void *u)"
void delete_hash_entry(void *p, void *q, void *u) {
  char *pp; TextHashTableVal *qq;
  pp= (char *)p;
  delete[] pp;
  qq= (TextHashTableVal *)q;
  delete qq;
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Adaptor for the Glib provided string hash function.
    (#g_str_hash#).
    @author M. Storti
    @param v (input) the value 
    @return the hash value provided for Glib
*/ 
unsigned int text_hash_func (const void *v) {
  const TextHashTableVal *vv = (TextHashTableVal *) v;
  return g_str_hash(vv->s);
}
#endif


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::TextHashTable ()"
TextHashTable::TextHashTable () {
  hash = g_hash_table_new(&g_str_hash,&g_str_equal);
};


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void TextHashTable::add_entry(char * key,char * value)"
void TextHashTable::add_entry(char * key,char * value) {
  TextHashTableVal *vold,*vnew;
  int keylen;
  keylen=strlen(key);
  vnew = new TextHashTableVal(value);
  vold = (TextHashTableVal *)g_hash_table_lookup(hash,key);
  if (vold) {
    printf("warning: redefining entry\n"
	   "key: %s\n"
	   "old value: %s\n"
	   "new value: %s\n",key,vold->s,vnew->s);
    delete vold;
  }
  // copy strings on new 
  char *keycp;
  keycp = local_copy(key);
  g_hash_table_insert (hash,keycp,vnew);
}

#if 0

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::get_entry_recursive(char * ,char *&, int&)"
void TextHashTable::get_entry_recursive(const char * key,const char *& svalue,
					int &glob_was_visited) {

  svalue = NULL;
  if (this==global_options) glob_was_visited=1;
  TextHashTableVal *value = (TextHashTableVal *) g_hash_table_lookup(hash,key);
  if (value) {
    if (value->called_times<INT_MAX) value->called_times++;
    svalue = value->s;
    return;
  }
  if (included_tables.size()==0) return;
  vector<TextHashTable *>::iterator k;
  for (k=included_tables.begin(); k!=included_tables.end(); k++) {
    (*k)->get_entry_recursive(key,svalue,glob_was_visited);
    if (svalue!=NULL) return;
  }
}
#endif


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::get_entry_recursive(const char *, const TextHashTableVal *&,int&)"
void TextHashTable::get_entry_recursive(const char * key,
					TextHashTableVal *& value,
					int &glob_was_visited) {
  if (this==global_options) glob_was_visited=1;
  value = (TextHashTableVal *)g_hash_table_lookup(hash,key);
  if (value!=NULL) return;
  if (included_tables.size()==0) return;
  vector<TextHashTable *>::iterator k;
  for (k=included_tables.begin(); k!=included_tables.end(); k++) {
    (*k)->get_entry_recursive(key,value,glob_was_visited);
    if (value!=NULL) return;
  }
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "" 
void TextHashTable::get_entry(const char * key,TextHashTableVal *& value) {
  int glob_was_visited=0;
  value=NULL;
  if (this) get_entry_recursive(key,value,glob_was_visited);
  //  if (value) return;
  if (!value && !glob_was_visited && global_options)
    global_options->get_entry_recursive(key,value,glob_was_visited);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "" 
void TextHashTable::get_entry(const char * key,const char *& svalue) {
  TextHashTableVal *value=NULL;
  get_entry(key,value);
  if (value && value->called_times<INT_MAX) value->called_times++;
  svalue = (value ? value->s : NULL);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::access_count()"
int TextHashTable::access_count(const char * key) {
  TextHashTableVal *value=NULL;
  get_entry(key,value);
  return (value ? value->called_times : -1);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::print(char * s = NULL ) const"
void TextHashTable::print(const char * s = NULL ) const {
  printf(((s == NULL) ? "Text hash table: \n" : "%s\n"),s);
  int count=0;
  g_hash_table_foreach (hash,&print_hash_entry,&count);
  if (count==0) printf("[No entries]\n");
  if (included_tables.size()==0) {
    printf("[No included hash tables]\n");
  } else {
    printf("Included hash tables:\n");
    for (int j=0; j < included_tables.size(); j++) {
      printf("%s  (%p)\n",included_tables_names[j]->c_str(),included_tables[j]);
    }
  }
  printf(" --\n");
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::~TextHashTable()"
TextHashTable::~TextHashTable() {
  g_hash_table_foreach (hash,&delete_hash_entry,NULL);
  g_hash_table_destroy (hash);
  for (int j=0; j<included_tables_names.size(); j++) {
    delete included_tables_names[j];
  }
  included_tables.resize(0);
  included_tables_names.resize(0);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::include_table(const string &s)"
void TextHashTable::include_table(const string &s) {
  THashTable::iterator k = thash_table.find(s);
  if (k==thash_table.end()) {
    printf("table not registered: %s\n",s.c_str());
    exit(1);
  }
  included_tables.push_back(k->second);
  string *sc= new string(s);
  included_tables_names.push_back(sc);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::register_name(const string &s)"
void TextHashTable::register_name(const string &s) {
  thash_table.insert(THashTableEntry(s,this));
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::read()"
void TextHashTable::read(FileStack *& fstack) {
  char *line;
  char *key, *val,*bsp=" \t";
  int he=0;
  while (1) {
    fstack->get_line(line);
    if (strstr("__END_HASH__",line)) break;
    key = strtok(line,bsp);
    val = strtok(NULL,"\n");
    if (!strcmp(key,"_table_include")) {
      include_table(val);
    } else {
      he++;
      if (val==NULL) val="";
      add_entry(key,val);
    }
  }
  g_hash_table_freeze(hash);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ ""
void TextHashTable::print_stat(void) {
  printf("------- TextHashTable statistics -------\n"
	 "[in brackets: number of times an entry was accessed]\n");
  string s;
  print_statistics=1;
  for (THashTable::iterator k=thash_table.begin();
       k!=thash_table.end(); k++) 
    {
      s = string("------ Table: \"");
      s.append(k->first.c_str()).append("\" -------");
      k->second->print(s.c_str());
    }
  print_statistics=0;
}
