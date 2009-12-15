#include "modulos.h"
#include <time.h>

void SlpTime(FILE *stream)

{
  static time_t t0=0;
  time_t t1;
  struct tm tm;

  t1=time(NULL);
  localtime_r(&t1,&tm);
  fprintf(stream,"Fecha actual: %d/%d/%d %d:%d:%d  ",
       tm.tm_mday,tm.tm_mon,tm.tm_year,tm.tm_hour,tm.tm_min,tm.tm_sec);
  if(t0 == 0) {
    t0=t1;
    fprintf(stream,"\n");
    }
  else {
    t1 -= t0;
    gmtime_r(&t1,&tm);
    fprintf(stream,"(Duracion: %d:%d:%d)\n",
      (tm.tm_mday-1)*24+tm.tm_hour,tm.tm_min,tm.tm_sec);
    t1 += t0;
    t0 = t1;
    }
  return;
}
