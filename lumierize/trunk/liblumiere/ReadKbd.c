#include "stdlib.h"
#include "readkbd.h"

#define  MAXSTR 1000

void kbdpause(void)
{
  char str[MAXSTR];
  setvbuf(stdin,"",_IOLBF,0);
  printf(" Press Enter to continue\n");
  getstrline(str,MAXSTR);
}

int readi(int n)
{
  char str[MAXSTR];
/*   //printf(" Asi %d\n",n); */
  printf("[%d]",n);
  setvbuf(stdin,"",_IOLBF,0);
  getstrline(str,MAXSTR);
  n= str[0]=='\0'?n:atoi(str);
  return(n);
}


double readd(double n)
{
  char str[MAXSTR];
  printf("[%11.7g]",n);
  setvbuf(stdin,"",_IOLBF,0);
  getstrline(str,MAXSTR);
/*   //printf(" Linea <<%s>>\n",str); */
  n= str[0]=='\0'?(double)n:(double)atof(str);
/*   printf(" h= %f  %f\n",n,atof(str)); */
  return(n);
}
float readf(float n)
{
  char str[MAXSTR];
  printf("[%9.5g]",n);
  setvbuf(stdin,"",_IOLBF,0);
  getstrline(str,MAXSTR);
  n= str[0]=='\0'?n:atof(str);
  return(n);
}
char readc (char c)
{
  char str[MAXSTR];
  printf("[%c]",c);
  setvbuf(stdin,"",_IOLBF,0);
  getstrline(str,MAXSTR);
  c= str[0]=='\0'?c:str[0];
  return(c);
}

void reads(const char* def,char *outstr)
{
   char str[MAXSTR]; 
/*  char *str; */
  printf("[%s] ",def);
/*   printf(" str %x\n",str); */
  setvbuf(stdin,"",_IOLBF,0);
  getstrline(str,MAXSTR);
  /*  str=readline(NULL); */
  /*   str[strlen(str)-1]='\0'; */
  /*   printf(" KAKA <%s>\n",str);   */
  /*    printf(" KAKA <%c>\n",*str);  */
  
/*   printf(" Sera aqui\n"); */
/*   printf(" str %x\n",str); */
  strcpy(outstr,(*str=='\0'?def:str));
/*   printf(" Esto es despues \n"); */
/*   printf(" str %x\n",str); */
/*   free(str); */
/*   printf(" DE aqui no paso\n"); */
}
 
