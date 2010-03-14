#include <stdlib.h>
#include <stdio.h>

#include "readkbd.h"

/* // Esta esta cogida de libro The C programming language . Por fin!! */
/* // El setvbuf lo que hace es limpiar el buffer del todo. Con eso */
/* // se consigue que no me quede ningun residuo que me haya dejado el  */
/* // scanf de los cojones. */
int getstrline(char s[], int lim)
{
  int c,i;
  setvbuf(stdin,"",_IOLBF,0);

  for(i=0;i<lim-1 && (c=getchar())!=EOF && c!='\n';++i) {
    s[i]=c;
  }
  s[i]='\0';
  return(i);
}

