#include "modulos.h"

static void makeurl_stsci(char *urldir, double ra, double dec, float xsize, float ysize);

int retrieveDSS(struct wcsimage *dssimage, double ra, double dec, float xsize, float ysize) {

  makeurl_stsci((*dssimage).image.file,ra,dec,xsize,ysize);
  printf(" URL %s\n",(*dssimage).image.file);
  ReadWCSImage(dssimage);
 

  return 0;

}

static void makeurl_stsci(char *urldir, double ra, double dec, float xsize, float ysize) {
  int rah, ram;
  float ras;
  int decg, decm;
  float decs;
  char dsig;

  r2hms(ra*M_PI/12, &rah, &ram, &ras);
  r2gms(dec*M_PI/180, &dsig,  &decg, &decm, &decs);

  sprintf(urldir,"http://archive.stsci.edu/cgi-bin/dss_search?v=1&RA=%02d%%20%02d%%20%05.2f&Dec=%%20%c%02d%%20%02d%%20%04.1f&e=J2000&h=%02d&w=%02d",rah,ram,ras,dsig,decg,decm,decs,(int)xsize,(int)ysize);

}
