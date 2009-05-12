#include "modulos.h"
#include <time.h>


void RNDInit(void)

{
  srand((unsigned int)time(NULL)/2);   /* // Inicio del generador random */
}
