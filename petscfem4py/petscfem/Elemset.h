// $Id: Elemset.h,v 1.1.2.8 2006/04/27 19:09:17 rodrigop Exp $

#ifndef PYPF_ELEMSET_H
#define PYPF_ELEMSET_H

#include <string>
#include <vector>
#include "petscfem4py.h"
#include "Object.h"


PYPF_NAMESPACE_BEGIN

class Elemset : SMARTPTR(Elemset)  
  public Object
{

 private:
  Elemset();
  
 protected:
  int nelem, nel;
  std::vector<int> icone;

 public:
  ~Elemset();
  Elemset(const Elemset& elemset);
  Elemset(const std::string& type,
	  int nelem, int nel, const int icone[]);

  std::string getType() const;
  std::string getName() const;
  void        setName(const std::string& name);

  void getData(int* nelem, int* nel, const int* icone[]) const;
  void setData(int  nelem, int  nel, const int  icone[]);
  void getDataSize(int* nelem, int* nel) const;

  void getPart(int* n, int* part[]) const;

  void getElem(int i, int* n, const int* elem[]) const;
  void setElem(int i, int  n, const int  elem[]);
  int  getSize() const;

  int  getNDof() const;
  void setNDof(int ndof);

public:
  void clear();
  void view() const;

};

PYPF_NAMESPACE_END

#endif // PYPF_ELEMSET_H

// Local Variables:
// mode: C++
// End:
