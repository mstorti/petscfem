/*__INSERT_LICENSE__*/
//$Id: idmap.cpp,v 1.3 2001/05/30 03:58:50 mstorti Exp $
 
#include <stdio.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <newmatio.h>
#include <cassert>
#include <math.h>

#include <sles.h>

#ifdef RH60
#include "libretto.h"
#else
#include <libretto/libretto.h>
#endif
#include <libretto/darray.h>

#undef HAVE_MEMMOVE // para que no chille al incluir petsccfonf.h
#include "texthash.h"
#include "fstack.h"
#include "utils.h"
#include "idmap.h"

using namespace std;

int n_created_rows=0,n_created_cols=0, n_deleted_rows=0, n_deleted_cols=0;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "axpy(row_t &y,double const a,const row_t &x)"
void axpy(row_t &y,double const a,const row_t &x) {
  row_t::const_iterator it,end_x,jy,end_y;
  end_x = x.end();
  end_y = y.end();
  for (it=x.begin(); it!= end_x; it++) {
    int j = it->first;
    double coef = it->second;
    double valy =0.;
    jy = y.find(j);
    if (jy != end_y) valy = jy->second;
    valy +=  a * coef;
    y[j] = valy;
  }
  erase_null(y);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "" 
void print(const row_t &y) {
  row_t::const_iterator it;
  for (it=y.begin(); it!=y.end(); it++) {
    printf("%d -> %f, \n",it->first,it->second);
  }
  printf("\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::get_block_matrix(int j,set<int> iindx,set<int> jindx,Matrix q)"
void idmap::get_block_matrix(int j,set<int> &iindx,set<int> &jindx,Matrix &qq) {

  vector<int> istack,jstack;
  set<int>::iterator it;
  row_t col,row;
  row_t::iterator itr,itc;

  jindx.insert(j);
  jstack.push_back(j);
  int ii,jj;

  while (jstack.size()>0) {

    // printf("loading row indices from col indices\n");

    // load all row indices connected to col indices in jstack
    while (jstack.size()>0) {
      jj = jstack.back();
      // printf("j index: %d\n",jj);
      jstack.pop_back();
      get_col(jj,col);
      for (itc=col.begin(); itc!=col.end(); itc++) {
	assert(itc->second != 0);
	ii = itc->first;
	if (iindx.find(ii) == iindx.end()) {
	  // printf("inserting i index: %d\n",ii);
	  istack.push_back(ii);
	  iindx.insert(ii);
	}
      }
    }

    // printf("loading col indices from row indices\n");

    // load all col indices connected to row indices in istack
    while (istack.size()>0) {
      ii = istack.back();
      // printf("row index: %d\n",ii);
      istack.pop_back();
      get_row(ii,row);
      for (itr=row.begin(); itr!=row.end(); itr++) {
	assert(itr->second != 0);
	jj = itr->first;
	if (jindx.find(jj) == jindx.end()) {
	  // printf("inserting col index: %d\n",jj);
	  jstack.push_back(jj);
	  jindx.insert(jj);
	}
      }
    }
  }    

  // load q matrix
  if (iindx.size()>0 && jindx.size()>0) {

    set<int>::iterator iti,itj;
    int mm,nn;
    mm=iindx.size();
    nn=jindx.size();

    qq.ReSize(mm,nn);
    int rindx=0;
    for (iti=iindx.begin(); iti!=iindx.end(); iti++) {
      rindx++;
      ii = *iti;
      int cindx=0;
      for (itj=jindx.begin(); itj!=jindx.end(); itj++) {
	cindx++;
	jj = *itj;
	double val;
	get_val(ii,jj,val);
	qq(rindx,cindx)= val;
      }
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::solve(double *x,double *y)" 
void idmap::solve(double *x,double *y) {

  int mm,nn;
  set<int> iindx,jindx;
  set<int>::iterator it;
  Matrix qq;
  ColumnVector xr,yr;
  set<int> done;

  for (int j=1; j<=n; j++) {
    VOID_IT(iindx);
    VOID_IT(jindx);
    if (find(done.begin(),done.end(),j)==done.end()) {
      get_block_matrix(j,iindx,jindx,qq);
      mm = iindx.size();
      nn = jindx.size();
      assert(mm>0 && nn>0);
      yr.ReSize(mm);
      xr.ReSize(nn);

      int ll=0;
      for (it=iindx.begin(); it!=iindx.end(); it++) {
	yr(++ll) = y[*it-1];
      }

      xr = (qq.t()*qq).i() * (qq.t()*yr);

//       cout << "qq: \n" << qq << endl;
//       cout << "yr: \n" << yr << endl;
//       cout << "xr: \n" << xr << endl;

      ll=0;
      for (it=jindx.begin(); it!=jindx.end(); it++) {
	  x[*it-1] = xr(++ll);
	  if (nn>1) done.insert(*it);
      }
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::row_set(const int i,const row_t &row)"
void idmap::row_set(const int i,const row_t &row) {

  del_row(i);
  row_t::const_iterator it;
  for (it=row.begin(); it!=row.end(); it++) {
    set_elem(i,it->first,it->second);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::del_col(const int j)"
void idmap::del_col(const int j) {

  row_t col;
  row_t::iterator it;
  get_col(j,col);
  for (it=col.begin(); it!=col.end(); it++) {
    set_elem(it->first,j,0.);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::column_set(const int j,const row_t &col)"
void idmap::column_set(const int j,row_t &col) {

  del_col(j);
  row_t::iterator it;
  for (it=col.begin(); it!=col.end(); it++) {
    set_elem(it->first,j,it->second);
  }
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::del_row(const int i)"
void idmap::del_row(const int i) {

  row_t row;
  row_t::iterator it;
  get_row(i,row);
  for (it=row.begin(); it!=row.end(); it++) {
    set_elem(i,it->first,0.);
  }
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::del_col(const int j)"
void idmap::del_col(const int j) {

  row_t col;
  row_t::iterator it;
  get_col(j,col);
  for (it=col.begin(); it!=col.end(); it++) {
    set_elem(it->first,j,0.);
  }
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::check()"
void idmap::check() {
  printf("checking consistency by rows: \n");
  row_t row,col;
  col_t *colp;
  row_t::iterator it;
  int j,i;

  int OK=1;

  for (i=1; i<=m; i++) {
    get_row(i,row);
    for (it=row.begin(); it!=row.end(); it++) {
      if (it->second==0.) {
	printf("(%d,%d) has zero value!!\n",i,it->first);
	OK=0;
      } else {
	j=it->first;
	get_col(j,col);
	if (!(col[i] == row[j])) {
	  printf("(%d,%d) -> %f has no column entry!!\n",i,j,it->second); 
	  OK=0;
	}
      }
    }
  }
#define is_OK (OK? "YES" : "NO")
  printf(" ------------> OK? : %s\n",is_OK);

  OK=1;
  printf("checking consistency by cols: \n");
  for (j=1; j<=n; j++) {
    get_col(j,col);
    for (it=col.begin(); it!=col.end(); it++) {
      i=it->first;
      get_row(i,row);
      if(row[j]==0) {
	printf("(%d,%d) -> %f has no row entry!!\n",i,j,it->second); 
	OK=0;
      }
    }
  }
  printf(" ------------> OK? : %s\n",is_OK);
#undef is_OK
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "void print(const col_t &row,char *s==NULL)" 
void print(col_t &col,char *s=NULL) {
  printf("%s", (s==NULL ? " " : s));
  col_t::iterator it;
  for (it=col.begin(); it!=col.end(); it++) {
    printf("%d ",*it);
  }
  printf("\n");
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void print(const row_t &row,char *s==NULL)" 
void print(row_t &row,char *s=NULL) {
  printf("%s", (s==NULL ? "" : s));
  row_t::iterator it;
  for (it=row.begin(); it!=row.end(); it++) {
    printf("(%d -> %f)  ",it->first,it->second);
  }
  printf("\n");
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "idmap::idmap(int m,map_type mtyp)" 
idmap::idmap(const int m_,map_type mtyp) {
  m=m_; n=m;
  ident = new int[m];
  iident = new int[n];
  row_map = new map<int,row_t*>;
  col_map = new map<int,col_t*>;

  if (mtyp == IDENTITY_MAP) {
    for (int k=1; k<=m; k++) {
      ident[k-1] = k;
      iident[k-1] = k;
    }
  } else {
    for (int k=1; k<=m; k++) {
      ident[k-1] = 0;
      iident[k-1] = 0;
    }
  }
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::get_val(const int i,const int j,double val)"
void idmap::get_val(const int i,const int j,double &val) {
  
  int jj=ident[i-1];
  if (jj==0) {
    val=0.;
  } else if (jj>0) {
    if (jj==j) {
      val=1.;
    } else {
      val=0.;
    }
  } else {
    row_t *row;
    row_t::iterator it;
    row = (*row_map)[i];
    it = row->find(j);
    if(it != row->end()) {
      val = (*row)[j];
    } else {
      val = 0.;
    }
  }
  return;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::print(char *s==NULL)"
void idmap::print(char *s=NULL) {

  row_t row;
  printf("%s",(s==NULL ? "" : s));
  for (int i=1; i<=m; i++) {
    get_row(i,row);
    if (row.size()>0) {
      printf("row %d: ",i);
      ::print(row);
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::print_by_col(char *s==NULL)"
void idmap::print_by_col(char *s=NULL) {

  printf("Size of matrix m,n: %d %d\n",m,n);
  row_t col;
  printf("%s",(s==NULL ? "" : s));
  for (int j=1; j<=n; j++) {
    get_col(j,col);
    if (col.size()>0) {
      printf("col %d: ",j);
      ::print(col);
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void erase_null(row_t &row)"
void erase_null(row_t &row) {
  row_t::iterator it;
  for (it=row.begin(); it!=row.end(); it++) {
    if (it->second==0.) row.erase(it);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void idmap::set_elem(const int i,const int j,const double val=0.)"
void idmap::set_elem(const int i,const int j,const double val=0.) {
  row_t row,col;
  col_t coll;
  get_row(i,row);
  get_col(j,col);
  if (val==0.) {
    row.erase(j);
    col.erase(i);
  } else {
    row[j] = val;
    col[i] = val;
  }

  // set row to NULL row
  if (ident[i-1]<0) {
    // row_map[i]->~row_t();
    delete (*row_map)[i];
    row_map->erase(i);
    n_deleted_rows++;
  }
    
  ident[i-1]=0;

  // clean row
  erase_null(row);

  // insert appropriate row
  if(row.size()==0) goto DONE;

  if (row.size()==1) {
    row_t::iterator it=row.begin();
    if (it->second==1.) {
      ident[i-1]=it->first;
      goto DONE;
    } 
  } 
  
  ident[i-1] = -1;
  (*row_map) [i] = new row_t;
  *(*row_map) [i] = row;
  n_created_rows++;
      
 DONE:   ;
  
  // set col to NULL
  if (iident[j-1]<0) {
    // assert(col_map->find(j) != col_map->end());
    delete (*col_map)[j];
    (*col_map).erase(j);
    n_deleted_cols++;
  }

  iident[j-1]=0;

  // clean col
  erase_null(col);

  // insert appropriate row
  if(col.size()==0) return;

  row_t::iterator it;
  if (col.size()==1) {
    it=col.begin();
    if (it->second==1.) {
      iident[j-1]=it->first;
      return;
    } 
  } 

  n_created_cols++;
  iident[j-1] = -1;
  // assert(col_map->find(j) == col_map->end());
#if 1
  (*col_map) [j] = new col_t;
#else
  pair<int,col_t*> pp;
  pp.first = j;
  pp.second = new col_t;
  col_map->insert(pp);
#endif
  for (it=col.begin(); it!=col.end(); it++) {
    (*col_map)[j] -> insert(it->first);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "idmap::get_col(const int j,row_t &col) " 
void idmap::get_col(const int j,row_t &col) {
  col_t *coll;
  int i = iident[j-1];
  VOID_IT(col);
  if (i>0) {

    col[i]=1.;

  } else if (i<0) {

    double val;
    int ii;
    coll = (*col_map)[j];
    col_t::iterator it;
    for (it=coll->begin(); it!=coll->end(); it++) {
      ii = *it;
      get_val(ii,j,val);
      col[ii] = val;
    }
  }
}
    
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "idmap::get_row(const int i,row_t &row) " 
void idmap::get_row(const int i,row_t &row) {
  int j = ident[i-1];
  VOID_IT(row);
  if (j>0) {
    row[j]=1.;
    return;
  } else if (j==0) {
    return;
  } else {
    row = *(*row_map)[i];
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "idmap::get_row(const int i,IdMapRow &row) " 
void idmap::get_row(const int i,IdMapRow &row) {
  int j = ident[i-1];
  if (j>0) {
    row.resize(1);
    row[0].j=j;
    row[0].coef=1.;
    return;
  } else if (j==0) {
    row.resize(0);
    return;
  } else {
    row_t::iterator it;
    row_t *row_;
    row_ = (*row_map)[i];
    row.resize(row_->size());
    int jj=0;
    for (it=row_->begin(); it!=row_->end(); it++) {
      row[jj].j = it->first;
      row[jj].coef = it->second;
      jj++;
    }
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//  Informe final sobre el remap_cols() bug:
//  ========================================
//  %
//  %
//  El efecto de este bug es que, despues de definir el idmap dofmap->id,
//  en readmesh, cuando hace la rutina `remap_cols()' que remapea las
//  columnas (grados de libertad) de forma de dejar las que corresponden
//  al procesador 0 primero, despues las del 1, etc... hasta nptoc-1, y
//  despues los fijos, entonces despues de remapear las columnas se morfa
//  una gran cantidad de memoria que despues no es liberada. No he llegado
//  a la raiz del problema pero he logrado los siguientes diagnosticos:
//  %
//  *./ La memoria perdida depende del tamanho maximo al cual ha llegado
//  los col_map y row_map y la suma de los row_t y col_t individuales, es
//  decir de todos los `map'. Esto parece ser un `bug' o por lo menos un
//  funcionamiento no deseado en el compilador o las librerias STL. 
//  %
//  *./ La memoria no se recupera ni siquiera destruyendo los objetos. Sin
//  embargo, si se vuelven a usar maps (eventualmente otros) la memoria es
//  reusada. Pareceria ser entonces que la implementacion de los maps
//  reserva una cantidad de memoria para para los maps que nunca es
//  liberada ni reducida, incluso si los maps son destruidos. 
//  %
//  *./ Hay un ejemplo tryme2.cpp que muestra esto.
//  %
//  *./ La solucion actual es cambiar el remap_cols() de manera que se van
//  cargando las filas del viejo idmap en un nuevo idmap y despues se
//  destruye el viejo. De esa forma la memoria perdida es nula porque no
//  se usan los maps ya que no hay columnas ni filas especiales. 
//  %
//  *./ Otra posibilidad seria reimplementar el idmap de otra forma, por
//  ejemplo tomando un `binary tree' de libretto, o algo mas pedestre
//  implementado por mi directamente. 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void idmap::remap_cols(int *perm) {
  row_t row;
  row_t::iterator l;

  // remaps rows
  for (int k=1; k<=m; k++) {
    get_row(k,row);
    del_row(k);
    for (l=row.begin(); l!=row.end(); l++) 
      set_elem(k,perm[l->first-1],l->second);
  }

#if 0
  wait_from_console("antes de copiar los row_map");  
  map<int,row_t*> *row_map_new = new map<int,row_t*>;
  *row_map_new = *row_map;

  map<int,col_t*> *col_map_new  = new map<int,col_t*>;
  *col_map_new = *col_map;
  wait_from_console("despues de copiar los row_map");  

  wait_from_console("antes de borrar los row_map");  
  delete row_map;
  delete col_map;
  wait_from_console("despues de borrar los row_map");  

  row_map = row_map_new;
  col_map = col_map_new;

#endif

  int n_new=0;
  for (int k=0; k<n; k++) 
    if (perm[k]>n_new)  n_new = perm[k];
  n = n_new;

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
idmap::~idmap(void) {

  row_map_t::iterator it;
  for (it=row_map->begin(); it!=row_map->end(); it++) 
    delete it->second;

  col_map_t::iterator itt;
  for (itt=col_map->begin(); itt!=col_map->end(); itt++) 
    delete itt->second;

  delete[] ident;
  delete[] iident;
  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void idmap::remap_cols(int *perm,idmap &idnew) {
  row_t row;
  row_t::iterator l;

  // remaps rows
  for (int k=1; k<=m; k++) {
    get_row(k,row);
    for (l=row.begin(); l!=row.end(); l++) 
      idnew.set_elem(k,perm[l->first-1],l->second);
  }

  int n_new=0;
  for (int k=0; k<n; k++) 
    if (perm[k]>n_new)  n_new = perm[k];
  idnew.n = n_new;

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void idmap::print_statistics(void) const{

  int nns=0,n_null=0;
  
  for (int k=0; k<m; k++) {
    if (ident[k]==-1) nns++;
    if (ident[k]==0) n_null++;
  }
  PetscPrintf(PETSC_COMM_WORLD,
	      "Number of null rows:       %d\n"
	      "Number of non-simple rows: %d\n"
	      "size of row_map:           %d\n",
	      n_null,nns,row_map->size());

  nns=0; n_null=0;
  for (int k=0; k<n; k++) {
    if (iident[k]==-1) nns++;
    if (iident[k]==0) n_null++;
  }
  PetscPrintf(PETSC_COMM_WORLD,
	      "Number of null cols:       %d\n"
	      "Number of non-simple cols: %d\n"
	      "size of col_map:           %d\n",
	      n_null,nns,col_map->size());

  PetscPrintf(PETSC_COMM_WORLD,
	      "Number of created rows,cols: %d, %d\n"
	      "Number of deleted rows,cols: %d, %d\n",
	      n_created_rows, n_created_cols,
	      n_deleted_rows, n_deleted_cols);
}
