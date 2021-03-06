//__INSERT_LICENSE__
//$Id: texthash.cpp,v 1.24.10.1 2007/02/19 20:23:56 mstorti Exp $
 
#include <iostream>
#include <sstream>
#include <string>
//  #include <string.h>
//  #include "util2.h"
#include <src/fem.h>
#include <src/texthash.h>

// The static member(s)
THashTable TextHashTable::thash_table;
TextHashTable *TextHashTable::global_options=NULL;

int TextHashTable::print_statistics=0;
TextHashTable *GLOBAL_OPTIONS;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTableVal::TextHashTableVal(char *)"
TextHashTableVal::TextHashTableVal(const char *s_) :
  s(local_copy(s_)), called_times(0) {}

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
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void TextHashTable::set_entries"
void TextHashTable
::set_entries(const std::map<std::string,std::string>& M) {
  using namespace std;
  map<string,string>::const_iterator m = M.begin();
  while (m != M.end()) {
    const char* key = m->first.c_str();
    const char* val = m->second.c_str(); 
    this->set_entry(key,val);
    m++;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void TextHashTable::set_entry"
void TextHashTable::set_entry(const char * key,const char * value) {
  TextHashTableVal *vold,*vnew;
  char *orig_key, *keycp;
  void *orig_key_v, *vold_v;
  int exists;
  // keylen=strlen(key);
  exists = g_hash_table_lookup_extended (hash,key,&orig_key_v,&vold_v);
  // vold = (TextHashTableVal *)g_hash_table_lookup(hash,key);
  if (exists) {
    vold = (TextHashTableVal *) vold_v;
    orig_key = (char *) orig_key_v;
    if (!strcmp(vold->s,value)) return;
    delete vold;
    keycp = orig_key; // reuse old key
  } else keycp = local_copy(key);

  vnew = new TextHashTableVal(value);
  // copy strings on new 
  g_hash_table_insert (hash,keycp,vnew);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void TextHashTable::add_entry"
void TextHashTable
::add_entry(const char * key,
	    const char * value,int warn) {
  TextHashTableVal *vold,*vnew;
  char *orig_key, *keycp;
  void *orig_key_v, *vold_v;
  int exists;
  // keylen=strlen(key);
  exists = g_hash_table_lookup_extended (hash,key,&orig_key_v,&vold_v);
  // vold = (TextHashTableVal *)g_hash_table_lookup(hash,key);
  if (exists) {
    vold = (TextHashTableVal *) vold_v;
    orig_key = (char *) orig_key_v;
    if (!strcmp(vold->s,value)) return;
    if (warn) printf("warning: redefining entry\n"
		     "key: %s\n"
		     "old value: %s\n"
		     "new value: %s\n",key,vold->s,value);
    delete vold;
    keycp = orig_key; // reuse old key
  } else keycp = local_copy(key);

  vnew = new TextHashTableVal(value);
  // copy strings on new 
  g_hash_table_insert (hash,keycp,vnew);
}

void TextHashTable
::add_entry(const char *key,const int *value,int n) {
  ostringstream s;
  for (int j=0; j<n; j++) s << value[j] << " ";
  s << ends;
  const char *val = s.str().c_str();
  add_entry(key,val);
  // s.freeze(0);
}

void TextHashTable
::add_entry(const char *key,const double *value,int n) {
  ostringstream s;
  for (int j=0; j<n; j++) s << value[j] << " ";
  s << ends;
  const char *val = s.str().c_str();
  add_entry(key,val);
  // s.freeze(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::get_entry_recursive"
void TextHashTable::get_entry_recursive(const char * key,
					TextHashTableVal *& value,
					int &glob_was_visited) const {
  if (this==global_options) glob_was_visited=1;
  // fixme:= cast to `(char *)' is for avoiding a warning with
  // old compiler versions
  value = (TextHashTableVal *)g_hash_table_lookup(hash,(char *)key);
  if (value!=NULL) return;
  if (included_tables.size()==0) return;
  vector<const TextHashTable *>::const_iterator k;
  for (k=included_tables.begin(); k!=included_tables.end(); k++) {
    (*k)->get_entry_recursive(key,value,glob_was_visited);
    if (value!=NULL) return;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "" 
void TextHashTable::get_entry(const char * key,TextHashTableVal *& value) const {
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
void TextHashTable::get_entry(const char * key,const char *& svalue) const {
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
void TextHashTable::print(const char * s) const {
  printf(((s == NULL) ? "Text hash table: \n" : "%s\n"),s);
  int count=0;
  g_hash_table_foreach (hash,&print_hash_entry,&count);
  if (count==0) printf("[No entries]\n");
  if (included_tables.size()==0) {
    printf("[No included hash tables]\n");
  } else {
    printf("Included hash tables:\n");
    for (unsigned int j=0; j < included_tables.size(); j++) {
      printf("%s  (%p)\n",included_tables_names[j]->c_str(),(void*)included_tables[j]);
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
  for (unsigned int j=0; j<included_tables_names.size(); j++) {
    delete included_tables_names[j];
  }
  included_tables.resize(0);
  included_tables_names.resize(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::include_table(const string &s)"
void TextHashTable
::include_table(const string &s,const TextHashTable *t) {
  const TextHashTable *tt;
  THashTable::iterator k;
  if (t==NULL) {
    k = thash_table.find(s);
    if (k==thash_table.end()) {
      printf("table not registered: %s\n",s.c_str());
      exit(1);
    }
    tt = k->second;
  } else {
    tt = t;
  }
  string *sc= new string(s);
  included_tables.push_back(tt);
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
  const char *key, *val;
  const char *bsp=" \t";
  int he=0;
  while (!fstack->get_line(line)) {
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
  PETSCFEM_ASSERT(fstack->last_error()==FileStack::read_ok,
		  "Couldn't read correctly hash table at\n"
		  "%s:%d: at (or after) line: \"%s\"",fstack->file_name(),
		  fstack->line_number(),
		  fstack->line_read());
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void TextHashTable::get_entry(const char *name,vector<double> &v) {
  const char *value;
  static const char *bsp=" \t";
  char *token;
  int k;
  double val;

  get_entry(name,value);
  if (value==NULL) return;

  char *buf= new char[strlen(value)+1];
  strcpy(buf,value);
  int j=0;
  while(1) {
    token = strtok(j++==0 ? buf : NULL ,bsp);
    if (!token) return;
    k = sscanf(token,"%lf",&val);
    PETSCFEM_ASSERT(k==1,
		    "Table entry does not contain a double\n"
		    "key: \"%s\", value: \"%s\", token \"%s\"\n",name,value,token);
    v.push_back(val);
  }
  delete[] buf;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
const TextHashTable * TextHashTable::find(const string &name) {
  THashTable::iterator k;
  k = thash_table.find(name);
  const TextHashTable *table = NULL;
  if (k!=thash_table.end()) table = k->second;
  return table;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

static void 
fill_std_map(void *p, void *q, void *m) {
  typedef TextHashTableVal THTVal;
  typedef std::map<std::string,std::string> MapStr;
  THTVal* thtval = reinterpret_cast<THTVal*>(q);
  MapStr& mapstr = *(reinterpret_cast<MapStr*>(m));
  const char* key = reinterpret_cast<const char*>(p);
  const char* val = reinterpret_cast<const char*>(thtval->s);
  mapstr[key] = val;
}

void TextHashTable::
get_entries(std::map<std::string,std::string>& M) const {
  void* m = reinterpret_cast<void*>(&M);
  g_hash_table_foreach (this->hash, &fill_std_map, m);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int remove_hash_entry(void *p, void *q, void *u)"
static int remove_hash_entry(void *p, void *q, void *u) {
  char *pp; TextHashTableVal *qq;
  pp = (char*) p;
  delete[] pp;
  qq = (TextHashTableVal*) q;
  delete qq;
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
#undef __FUNC__
#define __FUNC__ "TextHashTable::del_entries()"
void TextHashTable::del_entries() {
  g_hash_table_foreach_remove (hash, &remove_hash_entry, NULL);
  included_tables.resize(0);
  included_tables_names.resize(0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static void export1(void *p, void *q, void *ptr) {
  char *pp; 
  TextHashTableVal *qq;
  map<string,string> *table_p = (map<string,string> *)ptr;
  pp= (char *)p;
  qq= (TextHashTableVal *)q;
  (*table_p)[pp] = qq->s;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void TextHashTable
::export_entries(map<string,string> &table) {
  g_hash_table_foreach(hash,&export1,&table);
}
