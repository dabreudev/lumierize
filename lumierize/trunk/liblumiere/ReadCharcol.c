
/*---------------------------------------------------------------------------
  
  Version 19-February-1997                                    file:readnumcol.f
  @ ceg 
  ---------------------------------------------------------------------------
  Lee una columna de una fichero en forma de tabla 
  Devuelve una matriz.  
  DA 0 si sale satisfactoriamente 
  Da 1 si sale mal 
  
  Para compilar: 
  cc -c   $s2/Proced/C/modulos/ReadNumcol.c 
*/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "readcol.h"

#define DEBUG 0
int ReadCharcol(char file[],int col,char *vector,int *lvec,int charlen,int *nlin)
{
  
  int i;
  int a;
  char c;
  char word[1000];
  char nul[5000];
  int nl;
  FILE *fp;
  nl=0;

  if((fp=fopen(file,"r")) == NULL) {
    printf(" Cannot open file %s\n",file);
    exit(1);
  }
  while ( (a=getc(fp)) != EOF) 
    if(a == 10) {
      nl++; 
    }
  *nlin=nl;
  fclose(fp);
  fp=fopen(file,"r");
  for (i=0;i<nl;i++) {
    fgets(nul,5000,fp);
    if(DEBUG) printf(" nul <%s>\n",nul);
    c=nul[0];
    if(c=='#') {
      lvec[i]=0;
      vector[i*charlen]='\0';   
      continue;
    }
    LeeWord(nul,col,word);
    strncpy(vector+i*charlen,word,charlen);
    vector[i*charlen+charlen-1]='\0';  /*Aniadido ahora */
/*     if(DEBUG) printf(" vec <%s>\n",vector[i*charlen]); */
    if(DEBUG) printf(" vec <%s>\n",word);
    if(DEBUG) printf(" vec <%s>\n",vector+i*charlen);
    lvec[i]=1;

  }
  return 0 ;
} 

