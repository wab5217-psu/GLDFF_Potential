#include <stdlib.h>
#include <stdio.h>

int main(int argc,char *argv[]){

  char *fname=NULL;
  fname=argv[argc-1];

  char pot_file[256];
  
  sprintf(pot_file,"%s.pot",fname);
  fprintf(stderr,"%s\n",pot_file);

}
