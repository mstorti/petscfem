//__INSERT_LICENSE__
//$Id: getprop.cpp,v 1.5 2001/04/02 21:21:57 mstorti Exp $
  
#include "fem.h"
#include "readmesh.h"
#include "getprop.h"

#define MAXPROP 100

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
void print_prop_hash_entry(void *p, void *q, void *u) {
  char *pp; props_hash_entry *qq;
  pp= (char *)p;
  qq= (props_hash_entry *)q;
  printf("%s -> %d, %d\n",pp,qq->width,qq->position);
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void load_props(double *propel,int *elprpsindx,int nprops,double
		*elemprops) {
  for (int iprop=0; iprop<nprops; iprop++) {
    int indx=elprpsindx[iprop];
    if(indx!=-1) propel[iprop] = elemprops[indx];
    //    printf("prop %d = %f\n",iprop,propel[iprop]);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "get_double"
int get_double(TextHashTable *thash,const char *name,
	       double *retval,int defval=0,int n=1) {
  if (n==0) return 0;
  const char *value;
  char *token;
  int k,ierr;

  thash->get_entry(name,value);
  if (value==NULL) {
    if (defval) {
//        for (k=0; k<n; k++) {
//  	retval[k]=defval[k];
//        }
      return 0;
    } else {
      return 1;
    }
  }

  char *buf= new char[strlen(value)+1];
  strcpy(buf,value);
  for (k=0; k<n; k++) {
    token = strtok(k==0 ? buf : NULL ," ");
    PETSCFEM_ERROR(!token,
		   "Table entry does not contain enough data\n"
		   "key: %s, value: %s\n",name,value);
    sscanf(token,"%lf",&(retval[k]));
    PETSCFEM_ERROR(!token,
		   "Table entry does not contain a double\n"
		   "key: %s, value: %s\n",name,value);
  }
  delete[] buf;
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "get_int"
int get_int(TextHashTable *thash,const char *name,
	       int *retval,int defval=0,int n=1) {
  if (n==0) return 0;
  const char *value; 
  char *token;
  int k;

  thash->get_entry(name,value);
  if (value==NULL) {
    if (defval) {
//        for (k=0; k<n; k++) {
//  	retval[k]=defval[k];
//        }
      return 0;
    } else {
      return 1;
    }
  }

  char *buf= new char[strlen(value)+1];
  strcpy(buf,value);
  for (k=0; k<n; k++) {
    token = strtok(k==0 ? buf : NULL ," ");
    if (token==NULL && n==1) {
      // Act as get_flag, returns 1 if no value is assigned
      retval[k]=1;
      return 0;
    } else if (token==NULL) {
      return 1;
    }
    sscanf(token,"%d",&(retval[k]));
  }
  delete[] buf;
  return 0;
}

int get_string_from_string(string &buf,string &ret) {
  static regex_t regex[2];
  regmatch_t matchptr[2];
  
  // Compiles regexps only once.
  static int was_compiled=0;
  if (!was_compiled) {
    was_compiled=1;
    int ierr = regcomp (&regex[0],"^[ \t]*\"\\([^\"]*\\)\"",0);
    assert(ierr==0);
    ierr = regcomp (&regex[1],"^[ \t]*\\([^ \t]*\\)",0);
    assert(ierr==0);
  }

  for (int k=0; k<2; k++) {
    int flag = regexec(&regex[k],buf.c_str(),2,matchptr,0);
    if (flag==0) {
      int so=matchptr[1].rm_so;
      int eo=matchptr[1].rm_eo;
      ret = buf.substr(so,eo-so);
      buf.erase(so,eo-so);
      return 1;
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "get_string"
int get_string(TextHashTable *thash,const char *name,
	       string &ret,int defval=0,int n=1) {
  if (n==0) return 0;
  const char *value_;
  int k;

  // First of all, void the local buffer `ret'
//    if (ret!=NULL) {
//      delete[] ret;
//      ret=NULL;
//    }

  thash->get_entry(name,value_);
  if (value_ == NULL) {
    if (defval) {
      return 0;
    } else {
      return 1;
    }
  }

  string value = string(value_);
  int ierr=get_string_from_string(value,ret);
  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "get_prop"

/** This function returns the integer `iprop' which is the position in
    the vector propel[] where the user will find this property at each
    pass on the element loop. This property may be taken from the
    elemset hash table or from the per-element property table. In the
    first case, the property is put in propel[]. In the second case,
    the corresponding entry in propel is set to 0, and elprpspindx[iprop]
    is set to the corresponding column index in elemprops[]. 
*/
int get_prop(int & iprop,GHashTable *props,TextHashTable *thash,
	     int *elprpsindx,double *propel,const char *name,int n) {
  
  int ierr=0;
  props_hash_entry *phe;

  // looks int the properties-per-element table
  phe = (props_hash_entry *)g_hash_table_lookup(props,(void *)name);
  if(phe!=NULL) {
    //    printf("entry phe is: %d %d\n",phe->width,phe->position);
    int w = phe->width;
    if (n!=w) ierr = 1; CHKERRA(ierr); 
    for (int k=0; k<w; k++) {
      elprpsindx[iprop]=(phe->position)+k;
      propel[iprop]=0.;
      //      printf("set elprpsindx[%d] = %d\n",iprop,elprpsindx[iprop]);
      iprop++;
    }
    return 0;
  }

  // looks in the elemset table
  double val;
  ierr = get_double(thash,name,&val,1,n); CHKERRA(ierr);
  for (int k=0; k<n; k++) {
    elprpsindx[iprop]=-1;
    propel[iprop]=val;
    //    printf("set elprpsindx[%d] = %d\n",iprop,elprpsindx[iprop]);
    iprop++;
  }
  return 0;
}


/*
# Local Variables: $
# mode: c++ $
# End: $
*/
