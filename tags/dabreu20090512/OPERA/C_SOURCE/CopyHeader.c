#include "modulos.h"


int main() {


  fitsfile *f1,*f2;
  char  file1[200]="",file2[200]="";
  int status=0;
  int na=1;


  
  printf (" File to extract header: ");
  reads(file1,file1);
  printf (" File to write header: ");
  reads(file2,file2);
/*   printf(" Fichero 1 <<%s>>\n",file1); */
  
  if(ffopen(&f1, file1, READONLY, &status)) fits_report_error(stderr,status);
/*   exit(1); */
  if(ffinit(&f2, file2, &status)) fits_report_error(stderr,status);
/*   exit(1); */
/*   printf(" Hasta aqui bien\n"); */
  ffcphd(f1,f2,&status);
  ffuky(f2, TINT,"NAXIS1",&na,"Lenght of axis 1",&status);
  ffuky(f2, TINT,"NAXIS2",&na,"Lenght of axis 2",&status);
  na=0;
  ffuky(f2, TINT,"NAXIS",&na,"Number of axis",&status);
/*   printf(" Aqui falla\n"); */
  fits_close_file(f1,&status);
/*   printf(" Auiqn 1\n"); */
  fits_close_file(f2,&status);
/*   printf(" YA AA n"); */

  return 0;

}


