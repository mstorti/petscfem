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
            7, "steps", "serverhost", "port", "options", "step", "state_file", "record",
            3, "output_array_list", "output_field_list", "step");
    }
}
