/*
  This file belongs to the PETSc - FEM package, a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

#include <fstack.h>
#include <petsc.h>

int main () {
    char *line1,*line2;
    FileStack file1("file1.dat");
    FileStack file2("file2.dat");
    // FileStack file3("file3.dat");

  // Reads alternately one line from file1 and file2.
    int j=0;
    while (!file1.get_line(line1) && !file2.get_line(line2)) {
      printf("line %d on file_stack 1: <%s>\n",++j,line1);
      file1.unread_line("replaced line");
      file1.get_line(line1);
      printf("replaced line %d on file_stack 1: <%s>\n",++j,line1);
      printf("line %d on file_stack 2: <%s>\n",++j,line2);
    }
    file1.close();
    file2.close();

    // Reads from file2 into file1
    file1.open("file1.dat");
    file2.open("file2.dat");
    while (!file2.get_line(line2)) file1.unread_line(line2);
  
    printf("\n\n\ncat file1 reversed at the start of file2:\n");
    j=0;
    while (!file1.get_line(line1)) 
      printf("line %d : <%s>\n",++j,line1);
}
