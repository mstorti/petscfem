/*__INSERT_LICENSE__*/
//$Id: tfstack.cpp,v 1.6 2003/02/12 19:32:14 mstorti Exp $

#include <src/fstack.h>
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
    
    j=0;
    file1.close();
    file1.open("filen1.dat");
    FILE * fout = fopen("fstack2.out.tmp","w");
    while (!file1.get_line(line1)) 
      fprintf(fout,"%d>   %s:%d   <%s>\n",
	      ++j,file1.file_name(),file1.line_number(),line1);
    fclose(fout);
    file1.close();

    file1.open("file4.dat");
    j = 0;
    while (!file1.get_line(line1)) 
      fprintf(fout,"%d>   %s:%d   <%s>\n",
	      ++j,file1.file_name(),file1.line_number(),line1);
    if (file1.last_error()!=FileStack::eof) 
      printf("Couldn't process correctly file4.dat\n");
    file1.close();
}
