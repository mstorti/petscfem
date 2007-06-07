// $Id$

#ifndef PF4PY_DTABLE_H
#define PF4PY_DTABLE_H

#include <utility>
#include <algorithm>
#include <vector>
#include "namespace.h"
#include "Error.h"

#include "Object.h"

PF4PY_NAMESPACE_BEGIN

template<typename T>
class DTable
  : public Object
{
protected:
  std::pair<int,int> shape;
  std::vector<T>     array;

public:
  ~DTable() { }
  DTable() 
    : shape(0,0), array(0)
  { };
  DTable(const DTable& t) 
    : shape(t.shape), array(t.array) 
  { }
  DTable(int rows, int cols)
    : shape(rows,cols), array(rows*cols)
  { }
  DTable(int rows, int cols, const T array[])
    : shape(rows,cols), array(rows*cols)
  { 
    if (array != 0) 
      std::copy(array, array + rows*cols, &this->array[0]);
  }

public:
  inline operator T*() const
  { return const_cast<T*>(&this->array[0]); }

public:
  inline
  int getSize() const
  { return this->shape.first*this->shape.second; }

  inline 
  const std::pair<int,int>& getShape() const
  { return this->shape; }

  inline 
  const std::vector<T>& getArray() const
  { return this->array; };

};

PF4PY_NAMESPACE_END

#endif // PF4PY_DTABLE_H

// Local Variables:
// mode: C++
// End:
