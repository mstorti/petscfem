/* Automatically generated - may need to edit! */

#include <dx/dx.h>
#include <dx/modflags.h>

#if defined(intelnt)
#include <windows.h>
#endif

extern Error DXAddModule (char *, ...);

#if defined(intelnt)
void FAR WINAPI DXEntry()
#else
void DXEntry()
#endif
{
    {
        extern Error m_ExtProgImport(Object *, Object *);
        DXAddModule("ExtProgImport", m_ExtProgImport, 0,
            1, "socket_name",
            1, "value");
    }
}
