//__INSERT_LICENSE__
// $Id: fm2temp.cpp,v 1.1 2003/02/26 22:12:18 mstorti Exp $

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2Tmp::FastMat2Tmp() : size(0), store(NULL) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FastMat2::~FastMat2Tmp() { clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2Tmp::sync() {
  size = store_v.size();
  store = store_v.begin();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

FastMat2 & FastMat2::operator()(int j) {
  if (j>=size) { 
    store_v.resize(j+1); 
    for (int k=size; k<j; k++) store_v[k] = NULL;
    sync();
  }
  if (!store[j]) store[j] = new FastMat2;
  return *store[j];
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void FastMat2::clear() {
  for (int j=0; j<size; j++) delete store[j];
  store_v.clear();
  sync();
}
