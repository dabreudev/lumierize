
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
  cc -c   $s2/Proced/C/modulos/ReadWCScol.c 
*/
#include  "modulos.h"
int ReadWCScol(char file[],int col,double *vector,int *lvec,int *nlin)
{
  
  int i,j;
  int a;
  char c;
  char  word1[200];
  char nul[1000];
  int nl,column;
  
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
    column=col;
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

/*     //    printf(" 6666666666666666666666666666666\n"); */

/*     //  printf("Patatas %d\n",i); */
    LeeWord(nul,column,word1);
/*     //    printf("Patatas 2  %d\n",i); */
/*     //  printf("word 1 <<%s>>\n",word1); */
/*     //column++; */
/*     //    printf("Patatas 3  %d\n",i); */
/*     //LeeWord(nul,column,word2); */
/*     //  printf("word 2 <<%s>>\n",word2); */
/*     //column++; */
/*     //LeeWord(nul,column,word3); */
/*     //     printf("Patatas 4  %d\n",i); */
   

/*     //        printf("word 3 <<%s>>\n",word3); */
/*     //    printf("Patatas 5  %d\n",i); */

    if(!strcmp(word1,"INDEF")) {
      vector[i]=0.;
      lvec[i]=0;
      
    }
    else
      {
/* 	//		printf("aqui en medio %d\n",i); */

/* 	//vector[i]=atof(word1); */
/* 	//	printf(" 111 %f\n",vector[i]); */
/* 	//vector[i]=atof(word2); */
/* 	//	printf(" 222 %f\n",vector[i]); */
/* 	//vector[i]=atof(word3); */
/* 	//	printf(" 333 %f\n",vector[i]); */
/* 	//	printf("aqui en mas medio\n"); */
	lvec[i]=1;
        vector[i]=str2dec(word1);
/* 	//vector[i]=atof(word1)+atof(word2)/60.+atof(word3)/3600.;	 */
/* 	//exit(1); */
      }
/*     //      printf("Y ahora que\n"); */
  }
  return 0 ;
} 

