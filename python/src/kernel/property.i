
// %define %_property(Class, type, name, doc, getter, setter..)
// %wrapper {
// %define_as(%mangle(Class)_size_get(_t), _t->getSize() )
// }
// %extend Class { const int size; }
// %enddef


%define %property(Class, type, attr, doc, getcode, setcode...)
#if #set != ""
%wrapper { 
%#define %mangle(Class)_##attr##_get(self)      getcode } 
%#define %mangle(Class)_##attr##_set(self, val) setcode 
}
%extend Class { %arg(type) attr; }
#else
%wrapper { 
%#define %mangle(Class)_##attr##_get(self) getcode 
}
%extend Class { const %arg(type) attr; }
#endif
%enddef
