/// Class Layout for PETSc-FEM-2

class CellSet {
private:
  FieldList fl;
  CellSetList csl;
public:
  class FieldIndx;
  class CellSetIndx;
  class MapFun {
    void init(CellSet &cs,Task t);
    void proc(CellSet &cs,Task t,int cell);
    void close(CellSet &cs,Task t);
  };
  // Set and get values for particular cells
  // this cells should live in this processor,
  // otherwise a scatter is needed
  set(int cell,FieldIndx f,Scalar val);
  get(int cell,FieldIndx f);
  // Number of cells
  int size();
  // get a FieldIndx for a given FieldName
  FieldIndx get(char * field_name);
  // Add a field
  FieldIndx add(char *field_name);
  // Delete a field
  void del(char *field_name);
  // Manipulate connectivities / Cellsets
  CellSetIndx add(char *name,CellSet &cs,int n=1);
  CellSetIndx get(CellSet &cs,int j=0);
  CellSetIndx get(char *name,int j=0);
  // Manipulate cell connectivities
  void add_connect(int cell,CellSetindex csi,int cell_n);
  // MapFuns
  void add_mapfun(Task t,MapFun &f);
  void do_task(Task t);
};


