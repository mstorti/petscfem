// $Id$

#ifndef PF4PY_PROPTABLE_H
#define PF4PY_PROPTABLE_H

#include <utility>
#include <string>
#include <vector>
#include <set>

#include "namespace.h"
#include "Error.h"

#include "Object.h"
#include "DTable.h"


PF4PY_NAMESPACE_BEGIN

#define FieldEntry std::pair< std::string, int >
#define FieldList  std::vector< FieldEntry >

template<typename T>
class PTable
  : public Object
{

protected:
  inline void check_fields() const {
    typedef std::set<std::string> set_str;
    typedef set_str::iterator     set_str_it;
    set_str names;
    for (std::size_t i=0; i<this->fields.size(); i++) {
      const std::string& name = this->fields[i].first;
      std::pair<set_str_it,bool> p = names.insert(name);
      if (!p.second) 
	throw Error("PTable: duplicated field: '" + name + "'");
    }
  }
  inline void check_table() const {
    int cols = 0;
    for (std::size_t i=0; i<this->fields.size(); i++) {
      const int fs = this->fields[i].second;
      if (fs <= 0 ) throw Error("PTable: invalid field size" );
      cols += this->fields[i].second;
    }
    if (cols != this->table.getShape().second)
      throw Error("PTable: invalid number of columns in table");
  }
  inline void check() const { 
    this->check_fields(); 
    this->check_table();
  }

protected:
  FieldList  fields;
  DTable<T>  table;
  
public:
  PTable()
    : fields(),
      table()
      
  { };
  ~PTable() { }
  PTable(const PTable& t)
    : fields(t.fields),
      table(t.table)
  { }
  
  PTable(const FieldList& fields, const DTable<T>& table)
    : fields(fields),
      table(table)
  { this->check(); }

public:

  inline
  int getSize() const
  { return this->table.getSize(); }

  inline 
  const std::pair<int,int>& getShape() const
  { return this->table.getShape(); }

  inline 
  const std::vector<T>& getArray() const
  { return this->table.getArray(); };

  inline
  const FieldList& getFields() const
  { return this->fields; }

public:
  inline operator T*() const { return this->table; }

};

#undef FieldEntry
#undef FieldList

PF4PY_NAMESPACE_END

#endif // PF4PY_PROPTABLE_H

// Local Variables:
// mode: C++
// End:
