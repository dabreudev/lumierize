#include "modulos.h"


void stdastr(char *ar,char *dec,char *stdformat)
{

strncpy(stdformat+9,ar,2);
strncpy(stdformat+12,ar+2,2);
strncpy(stdformat+15,ar+4,4);

strncpy(stdformat+22,dec,3);
strncpy(stdformat+26,dec+3,2);
strncpy(stdformat+29,dec+5,2);

}
