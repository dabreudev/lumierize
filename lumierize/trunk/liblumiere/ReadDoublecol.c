
/*---------------------------------------------------------------------------
  
  Version 19-February-1997                                    file:readnumcol.f
  @ ceg 
  ---------------------------------------------------------------------------
  Lee una columna de una fichero en forma de tabla 
  Devuelve una matriz. Las lineas que comiencen por # las ignora 
  Los valores que contengan INDEF los pone a 0. 
  DA 0 si sale satisfactoriamente 
  Da 1 si sale mal, por ejemplo si algunos de los valores no es un numero 
  
  Para compilar: 
  cc -c   $s2/Proced/C/modulos/ReadNumcol.c 
*/
#include <string.h>
#include <stdio.h>
#include "readcol.h"
int ReadDoublecol(char file[],int col,double *vector,int *lvec,int *nlin)
{
  
  int i,j;
  int a;
  char c;
  char  word[200];
  char nul[1000];
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
/*     //    printf("ReadNumcol reading %d row",i); */
    fgets(nul,1000,fp);
    c=nul[0];
    if(c=='#') {
      lvec[i]=0;
      continue;
    }
    j=0;
    while(j<strlen(nul)) {
      if(nul[j]=='\t') {
	nul[j]=' ';
      }
      j++;
    }
/*     //    printf("Pasa \n"); */
    strtok(nul,"\n");
    strcat(nul," ");
    LeeWord(nul,col,word);
/*     //    printf("word <<%s>>\n",word); */
    
    if(!strcmp(word,"INDEF")) {
      vector[i]=0.;
      lvec[i]=0;
      
    }
    else
      {
/* 	//printf("aqui en medio\n"); */
	vector[i]=(double)atof(word);
/* 	//printf("aqui en mas medio\n"); */
	lvec[i]=1;
      }
/*     //  printf("Y ahora que\n"); */
  }
  return 0 ;
} 

