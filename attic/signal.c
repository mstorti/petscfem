#define _GNU_SOURCE

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>

int stop=0;     

void action(int sig) { 
  stop++;
  if (stop>1) {
    signal(SIGINT,SIG_DFL);
    raise(SIGINT);
  }
}
     
int main (void) {
  int j=0;
  char ans;

  signal(SIGINT,&action);
  while(1) { 
    printf("j: %d\n",j++);
    if (stop) {
      printf("Exit ? [y/n] > ");
      fflush(stdout);
      scanf("%c",&ans);
      if (ans=='y') {
	printf("Exiting by user request...\n");
	exit(1);
      }
      while (ans!='\n') scanf("%c",&ans);
      stop=0;
    }
    sleep(1);
  }
}
