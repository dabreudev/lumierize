#include "modulos.h"

int f_kfl(FILE *fp,char *pclave, char label[][71])

{
  char bloq[80];
  int i,j=0;

  rewind(fp);
  f_clrl(label);

  do {
    i=fread(bloq,1,80,fp);
    if(strstr(bloq,pclave) != NULL) {
      if(!strcmp(pclave,"COMMENT ")) strncpy(label[j],bloq+10,70);
      if(!strcmp(pclave,"HISTORY ")) strncpy(label[j],bloq+10,70);
      j++;
      }
    } while ((j < MAX_KEYLAB)&&(strstr(bloq,"END      ") == NULL)&&(i != 0));
  rewind(fp);

  if(!j) printf("f_kfl: WARNING. No he encontrado la Keyword %s\n",pclave);
  if(j == MAX_KEYLAB) 
    printf("f_kfl: WARNING. No. maximo de comentarios encontrados (%d)\n",
           MAX_KEYLAB);
  return(j);
}
