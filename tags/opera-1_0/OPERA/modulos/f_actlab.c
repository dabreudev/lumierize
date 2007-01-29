#include "modulos.h"

int f_actlab(char label[][71], char *var)

{
  int i=0;
 
  if(strlen(var) > 70) {
    fprintf(stderr,"f_actlab: ERROR. El var debe tener menos de 70 caract\n");
    return(0);
    }
 
 
  for(i=0; i<MAX_KEYLAB; i++) if(label[i][0]==' ') break; 

  if(i==MAX_KEYLAB) {
    i=MAX_KEYLAB-1;
    memcpy(label[0],label[1],71*(MAX_KEYLAB-1));
    }

  strcpy(label[i],var);

  return(1);
 }
