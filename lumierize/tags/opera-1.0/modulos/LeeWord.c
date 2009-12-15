#include "modulos.h"

void LeeWord(char a[],int nw,char word[])
{
  char *b;
  char d[5000]; 
  int i=0;
  b=a;      
  
/*   printf(" Llamandao e a ae \n"); */

  while(i<nw) {
    while(strcspn(b," ")==0) b++; 
    strcpy(d,b);
    strtok(d," ");
    i++;
    b=strpbrk(b+1," ");
    
  }
/*       printf("La d >>%s<< %d\n",d, nw); */
    strcpy(word,d);
/*     memcpy(word,d,1000); */
/*     printf(">>%s<<\n",word); */

}



