//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

class Elemset;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

void bless_elemset0(char*, Elemset *&);
void bless_elemset_ns(char*, Elemset *&);

void bless_elemset(char *type,Elemset *& elemset) {
  elemset = 0;
  bless_elemset0(type,elemset);
  if (elemset) return;
  bless_elemset_ns(type, elemset);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

int ns_main(int,char**);

int main(int argc,char **args) 
{
  return ns_main(argc,args);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
