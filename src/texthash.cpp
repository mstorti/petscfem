/*
  This file belongs to the PETSc - FEM package, a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/
 
#include <stdio.h>
#include <string.h>
#include "texthash.h"

// The static member(s)
THashTable TextHashTable::thash_table;
TextHashTable *TextHashTable::global_options=NULL;


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void print_hash_entry(void *p, void *q, void *u)"
void print_hash_entry(void *p, void *q, void *count) {
  char *pp; char *qq;
  pp= (char *)p;
  qq= (char *)q;
  printf("%s -> %s\n",pp,qq);
  (*(int *)count)++;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void delete_hash_entry(void *p, void *q, void *u)"
void delete_hash_entry(void *p, void *q, void *u) {
  char *pp; char *qq;
  pp= (char *)p;
  qq= (char *)q;
  delete[] pp;
  delete[] qq;
}

#if 0

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "unsigned int hash_func(const void *kkey)"
unsigned int hash_func(const void *kkey) {
  unsigned int qq;
  size_t ilen = sizeof(int);
  qq = 0;
  char *key = (char *)kkey;
  char *cqq=(char *)(&qq);
  int len = strlen((char *)key);
  len = (strlen(key)>sizeof(int) ? sizeof(int) : len);
  for (int k=0; k<len; k++) {
    strncpy(cqq+len-k-1,key+k,1);
  }
  return qq;
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int key_compare_func(const void *key1, const void *key2)"
int key_compare_func(const void *key1, const void *key2) {
  char *q1,*q2;
  q1=(char *)key1;
  q2=(char *)key2;
  return !strcmp(q1,q2);
}
#endif


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::TextHashTable ()"
TextHashTable::TextHashTable () {
  hash = g_hash_table_new(&g_str_hash,&g_str_equal);
  // hash = g_hash_table_new(&hash_func,&key_compare_func);
};


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void TextHashTable::add_entry(char * key,char * value)"
void TextHashTable::add_entry(char * key,char * value) {
  int keylen,vallen;
  keylen=strlen(key);
  vallen=strlen(value);
  char *p = (char *)g_hash_table_lookup(hash,key);
  if (p) {
    printf("warning: redefining entry\n"
	   "key: %s\n"
	   "old value: %s\n"
	   "new value: %s\n",key,p,value);
    delete[] p;
  }
  // copy strings on new 
  char *keycp, *valcp;
  keycp = new char[keylen+1];
  valcp = new char[vallen+1];
  strcpy(keycp,key);
  strcpy(valcp,value);
  //  printf("insert %s -> \"%s\"\n",keycp,valcp);
  g_hash_table_insert (hash,keycp,valcp);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::get_entry(char * key,char *& value)"
void TextHashTable::get_entry_recursive(char * key,char *& value,
					int &glob_was_visited) {

  if (this==global_options) glob_was_visited=1;
  value = (char *)g_hash_table_lookup(hash,key);
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
void TextHashTable::get_entry(char * key,char *& value) {
  int glob_was_visited=0;
  value=NULL;
  if (this) get_entry_recursive(key,value,glob_was_visited);
  if (value) return;
  if (!glob_was_visited && global_options)
    global_options->get_entry_recursive(key,value,glob_was_visited);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "TextHashTable::print(char * s = NULL ) const"
void TextHashTable::print(char * s = NULL ) const {
  printf(((s == NULL) ? "Text hash table: " : "%s\n"),s);
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
  void * dummy;
  printf("deleting TextHashTable!!!\n");
  g_hash_table_foreach (hash,&delete_hash_entry,dummy);
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
