#include <stdio.h>
#include "filenlin.h"


int FileNLin(char *file)

{
  FILE *fp;
  int n=0;
  int a;

  if( (fp = fopen(file,"r")) == NULL ) {
    fprintf(stderr,"FileNLin: ERROR. No puedo abrir el fichero %s\n",file);
    return(0);
    }

  while ( (a=getc(fp)) != EOF)
    if(a == 10) n++;
  fclose(fp);
  return(n);
}

