
#include "modulos.h"

/* // Esta esta cogida de libro The C programming language . Por fin!! */
/* // El setvbuf lo que hace es limpiar el buffer del todo. Con eso */
/* // se consigue que no me quede ningun residuo que me haya dejado el  */
/* // scanf de los cojones. */
int getline(s,lim)
     char s[];
     int lim;
{
  int c,i;
  setvbuf(stdin,"",_IOLBF,0);

  for(i=0;i<lim-1 && (c=getchar())!=EOF && c!='\n';++i) {
    s[i]=c;
  }
  s[i]='\0';
  return(i);
}


/* //Esta es igual pero hecha por mi. Es para leer de un fichero una linea */
int fgetline(stream, s, lim)
     FILE *stream;
     char s[];
     int lim;
{
  int c,i;
  setvbuf(stream,"",_IOLBF,0);

  for(i=0;i< (c=fgetc(stream))!=EOF && c!='\n';++i) {
    s[i]=c;
  }
  s[i]='\0';
  return(i);
}

