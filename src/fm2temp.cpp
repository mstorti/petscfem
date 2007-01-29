//__INSERT_LICENSE__
// $Id: fm2temp.cpp,v 1.3.92.1 2007/01/29 21:07:57 dalcinl Exp $

#include <src/fastmat2.h>
#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2Tmp::FastMat2Tmp() : store(NULL) , size(0) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2Tmp::~FastMat2Tmp() { clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2Tmp::sync() {
  size = store_v.size();
  store = &*store_v.begin();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

FastMat2 & FastMat2Tmp::operator()(int j) {
  if (j>=size) { 
    store_v.resize(j+1); 
    for (int k=size; k<j; k++) store_v[k] = NULL;
    sync();
  }
  if (!store[j]) store[j] = new FastMat2;
  return *store[j];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2Tmp::clear() {
  for (int j=0; j<size; j++) delete store[j];
  store_v.clear();
  sync();
}
