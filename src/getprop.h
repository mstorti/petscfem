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
 

#ifndef GETPROP_H
#define GETPROP_H

#include <regex.h>
#include <string>

#include "texthash.h"

/**@name getprop package */
//@{

#define DEFPROP(name) \
  int name##_indx = iprop; \
  ierr = get_prop(iprop,elem_prop_names,thash,elprpsindx,propel, \
		  #name,1); \

#define USE_DEFAULT 1
#define DONT_USE_DEFAULT 0

/** Fast load of element properties.  It transparently loads
    properties defined `per element' or globally for the elemenset in
    the properties hash table.  You should call DEFPROP(prop1);
    DEFPROP(prop2); etc.. before entering the element loop, and inside
    the element loop you call load\_props(...), then you can use the
    values as propel[prop1\_indx], propel[prop2\_indx], etc...

    @author M. Storti
    @param propel a double working array defined by the user.
    @param elprpsindx integer working array defined by the user.
    @param nprops number of properties to be `fast loaded'
    @param elemprops array of `per element' properties
 */
void load_props(double *propel,int *elprpsindx,int nprops,double *elemprops);

/** Gets a double value from the hash table.

    @author M. Storti
    @param thash pointer to the elemset hash table.
    @param name name of the variable
    @param retval the value read from the table. Set to the default
    value if `defval=0'.
    @param n number of doubles to be read
*/
int get_double(TextHashTable *thash,const char *name,
	       double *retval,int defval=0,int n=1);

// fixme:= REWRITE!!!!!!!!!!!!
/** Gets a string value from the hash table.

    @author M. Storti
    @param thash pointer to the elemset hash table.
    @param name name of the variable
    @param retval the value read from the table. Set to the default
    value if `defval=0'. Strings may be enclosed in double quotes to further
    @param n number of integers to be read
*/
int get_string(TextHashTable *thash,const char *name,
	       string &ret,int defval=0,int n=1);

/** Gets an integer value from the hash table.
    @author M. Storti
    @param thash pointer to the elemset hash table.
    @param name name of the variable
    @param retval the value read from the table. Set to the default
    value if `defval=0'. If the entry is found but not value assigned,
    then assign 1 (to be used as a `get\_flag()' function). 
    @param n number of integers to be read
*/
int get_int(TextHashTable *thash,const char *name,
	       int *retval,int defval=0,int n=1);

#if 0
int get_flag(TextHashTable *thash,const char *name,
	       int *retval,int defval=0,int n=1);
#endif

/** Prepares for a `fast load' of properties. 
    @author M. Storti
*/
int get_prop(int & iprop,GHashTable *props,TextHashTable *thash,
	     int *elprpsindx,double *propel,const char *name,int n);
//@}

#endif
