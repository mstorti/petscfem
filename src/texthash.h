// -*- mode: c++ -*-

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

#ifndef __TEXTHASH_H__
#define __TEXTHASH_H__

#include <map>
#include <vector>
#include <string>
#include <cstring>

#include <glib.h>

#include "fstack.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Makes a temporary copy of a string.
    @author M. Storti
    @param cstr (input) the string to be copied
    @return a pointer to the copied string
*/ 
char * local_copy(const char * cstr);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** @name Text hash functions */
//@{
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** To be passed to Glib traversal functions to print the entire
    hash. 
    @author M. Storti
    @param p (input) key string
    @param q (input) value string
    @param u (not used) as required by Glib
*/ 
void print_hash_entry(void *p, void *q, void *u);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Remove entry from hash. To be passed to Glib
    traversal functions to build the hash destructor. 
    @author M. Storti
    @param p (input) key string
    @param q (input) value string
    @param u (not used) as required by Glib
*/ 
void delete_hash_entry(void *p, void *q, void *u);

//@}

class TextHashTable;
typedef pair<string,TextHashTable *> THashTableEntry;
typedef multimap<string,TextHashTable *> THashTable;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Text hash tables (key and value are strings) are used to store
    elemset properties. 
*/ 
class TextHashTable {
public:
  /** Prints the entire hash. 
      @author M. Storti
      @param s (input) optional string
  */ 
  void print(const char * = NULL) const;

  /** Adds an entry to the hash.
      @author M. Storti
      @param key (input) key of the entry
      @param value (input) value of the entry
  */ 
  void add_entry(char *key,char *value);

  /** Adds an included table. 
      @author M. Storti
      @param ithash (input) the table to be included
      @param s (input) its name
  */ 
  void include_table(const string &s);

  /** Registers the table by its name.
      @author M. Storti
      @param s (input) the name of the table. 
  */
  void register_name(const string &s);

  /** Sets this table as the global one.
      @author M. Storti
  */
  void set_as_global() {global_options=this;};

  /** Searches an entry in the hash. 
      @author M. Storti
      @param key (input) key of the entry
      @param value (output) value of the entry
  */ 
  void get_entry(const char *,const char *&);

  /** Returns the number of times a particular key was accessed.
      @author M. Storti
      @param key (input) key of the entry
      @return the number of access to the key
  */ 
  int access_count(const char *);

  /** Constructs a void hash table.  
      @author M. Storti
  */ 
  TextHashTable();

  /// Destructor. 
  ~TextHashTable();

  /** Reads a text hash table from a filestack. 
      @author M. Storti
      @param fstack (input) The filestack from which the hash table is
      read. 
  */ 
  void read(FileStack *& fstack);

  /** Print all the text-hash-tables
      @author M. Storti
  */ 
  static void print_stat();

private:
  /// The underlying Glib hash. 
  GHashTable *hash;  

  /// A list of pointers to other (included) hashes
  vector<TextHashTable *> included_tables;

  /// A list of the names of the included hashes
  vector<string *> included_tables_names;

  /// The global register table
  static THashTable thash_table;

  /// The global hash table
  static TextHashTable *global_options;

  /** Searches an entry in the hash recursively. 
      Returns the whole entry (a struct) instead of the plain string.
      @author M. Storti
      @param key (input) key of the entry
      @param value (output) value of the entry
  */ 
  void get_entry_recursive(const char *,TextHashTableVal *&,
			   int &glob_was_visited);

  /** Searches an entry in the hash. 
      This returns the whole entry (a struct) instead 
      of the plain string.
      @author M. Storti
      @param key (input) key of the entry
      @param value (output) value of the entry
  */ 
  void get_entry(const char *,TextHashTableVal *&);
};

#endif /* __TEXTHASH_H__ */
