// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: texthash.h,v 1.7.4.1 2001/12/26 15:36:13 mstorti Exp $

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
// This is a typical object in the table of thash's
typedef pair<string,TextHashTable *> THashTableEntry;
// This is the global table of thash's
typedef multimap<string,const TextHashTable *> THashTable;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** TextHashTable's are a map string -> TextHashTableVal objects
*/ 
class TextHashTableVal {
public:
  // string containing the value for the option
  char *s;
  // Number of times this option was accessed for statistic purposes
  int called_times;
  // Constructor from the string
  TextHashTableVal(const char *s_=NULL);
  // Destructor
  ~TextHashTableVal() {delete s;};
};

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
  void add_entry(const char *key,const char *value);

  /** Adds an entry to the hash.
      @author M. Storti
      @param key (input) key of the entry
      @param value (input) value of the entry
  */ 
  void add_entry(const char *key,const int *value,int n=1);
  /** Adds an entry to the hash.
      @author M. Storti
      @param key (input) key of the entry
      @param value (input) value of the entry
  */ 
  void add_entry(const char *key,const double *value,int n=1);

  /** Adds an included table. 
      @author M. Storti
      @param s (input) its name
      @param t (input) a pointer to the TextHashTable to be included
      (optional). If not given, look for `s' in the
      `included_tables_names' table
  */ 
  void include_table(const string &s,const TextHashTable *t = NULL);

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

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** To be passed to Glib traversal functions to print the entire
      hash. 
      @author M. Storti
      @param p (input) key string
      @param q (input) value string
      @param u (not used) as required by Glib
  */ 
  friend void print_hash_entry(void *p, void *q, void *u);

private:
  /// The underlying Glib hash. 
  GHashTable *hash;  

  /// A list of pointers to other (included) hashes
  vector<const TextHashTable *> included_tables;

  /// A list of the names of the included hashes
  vector<const string *> included_tables_names;

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
			   int &glob_was_visited) const;

  /** Searches an entry in the hash. 
      This returns the whole entry (a struct) instead 
      of the plain string.
      @author M. Storti
      @param key (input) key of the entry
      @param value (output) value of the entry
  */ 
  void get_entry(const char *,TextHashTableVal *&);

  static int print_statistics;
};

#endif /* __TEXTHASH_H__ */
