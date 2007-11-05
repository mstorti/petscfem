//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

class Elemset;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

void bless_elemset0(char*, Elemset *&);
void bless_elemset_advdif(char*, Elemset *&);

void bless_elemset(char *type,Elemset *& elemset) {
  elemset = 0;
  bless_elemset0(type,elemset);
  if (elemset) return;
  bless_elemset_advdif(type, elemset);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

int advdif_main(int,char**);

int main(int argc,char **args) 
{
  return advdif_main(argc,args);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
