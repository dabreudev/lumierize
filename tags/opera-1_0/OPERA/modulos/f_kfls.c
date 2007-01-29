#include "modulos.h"

int f_kfls(FILE *fp,char *pclave, char *var)

{
  char *p,bloq[2880];
  int i,j;

  rewind(fp);
  do {
    i=fread(bloq,1,2880,fp);
    if((p=strstr(bloq,pclave)) != NULL) {
      if(*(p+10) == '\'')  {
	j=0;
	while(*(p+11+j) != '\'' && j < 65) j++;
        strncpy(var,p+11,j);
	var[j]='\0';
        }
      }
    } while ((p == NULL) && (strstr(bloq,"END      ") == NULL) && (i != 0));
  rewind(fp);

  if(p!=NULL) return(1);
  else {
    printf("f_kfs: WARNING. No he encontrado la Keyword %s\n",pclave);
    return(0);
    }
}
