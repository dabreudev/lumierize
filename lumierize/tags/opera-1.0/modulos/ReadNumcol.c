
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
#include  "modulos.h"
#include <sys/time.h>
#include <unistd.h>

#define DEBUG 0
#define DEBUG2 0


int ReadNumcol(char file[],int col,float *vector,int *lvec,int *nlin)
{
  
  int i,j;
  int a;
  char c;
  char  word[500];
  char nul[5000];
  int nl;

  int nullen;
  
  
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

    
    
    /*         printf("ReadNumcol reading %d row",i); */
    fgets(nul,5000,fp);

    
    if(DEBUG)  printf(" NUL <%s>\n",nul);
    c=nul[0];
    if(c=='#') {
      lvec[i]=0;
      continue;
    }
    j=0;
    if(DEBUG2) printf(" lon %d\n",strlen(nul));

    
    nullen=strlen(nul);

    while(j<nullen) {
      if(nul[j]=='\t') {
	nul[j]=' ';
      }
      j++;
    }
    

    
if(DEBUG)      printf("Pasa \n"); 
    strtok(nul,"\n");
    strcat(nul," ");
if(DEBUG)  printf(" NUL2 <%s>\n",nul);

    LeeWord(nul,col,word);


if(DEBUG)    printf("word <<%s>>\n",word); 
    
    if(!strcmp(word,"INDEF")) {
      vector[i]=0.;
      lvec[i]=0;
      
    }
    else
      {
/* 	//	printf("aqui en medio\n"); */
	vector[i]=atof(word);
/* 	//	printf("aqui en mas medio\n"); */
	lvec[i]=1;
      }
/*     //  printf("Y ahora que\n"); */



  }
  return 0 ;
} 

