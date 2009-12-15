#include "modulos.h"


struct calib_files {
  struct image dark1d;
  struct image flat1d;
  struct image errdark1d;
  struct image errflat1d;
  float gain;
  float ron;
  char trimsec[100];
};

struct science_files{

  struct image raw;
  struct image errraw;
  struct image reduced;
  struct image errreduced;
  struct image trimreduced;
  struct image errtrimreduced;
};

struct calib_files CAL;
struct science_files SCI;

void LoadParam_kbd(void);
void ReadInputImages(void);
void PerformOperations(void);
void InitOuputImages(void);
void CheckDim(void);
void TrimImages(void);
void SaveOutputImages(void);
void CloseFiles(void);


int main()
{
/*   LoadParam_kbd(); */
/*   TrimImages(); */
/*   exit(1); */
/*   int i=0; */

/*   cpgbeg(0,"?",1,1);   */
  LoadParam_kbd();
  ReadInputImages();

  CheckDim();
  InitOuputImages();
  PerformOperations();

  TrimImages();

  SaveOutputImages();
  CloseFiles();

  return(0);
 
}




void LoadParam_kbd(void) {

  printf("\n Input FITS file to reduce: ");
  reads("",SCI.raw.file);
  printf("\n Input FITS output reduced file: ");
  reads("",SCI.trimreduced.file);
  printf("\n Input FITS output error reduced file: ");
  reads("",SCI.errtrimreduced.file);
  printf("\n Input FITS file with 1-d DARK: ");
  reads("",CAL.dark1d.file);
  printf("\n Input FITS file with 1-d error DARK: ");
  reads("",CAL.errdark1d.file);
  printf("\n Input FITS file with 1-d FLAT: ");
  reads("",CAL.flat1d.file);
  printf("\n Input FITS file with 1-d error FLAT: ");
  reads("",CAL.errflat1d.file);
  printf(" Input TRIM section (Iraf style): ");
  reads("",CAL.trimsec);
  printf(" Input Gain: ");
  CAL.gain=readf(1.);
  printf(" Input RON: ");
  CAL.ron=readf(1.);
}


void ReadInputImages() {
  ReadImage(&SCI.raw);
  ReadImage(&CAL.dark1d);
  ReadImage(&CAL.errdark1d);
  ReadImage(&CAL.flat1d);
  ReadImage(&CAL.errflat1d);

  printf(" npix SCI.raw %ld \n",SCI.raw.npixels);
  printf(" npix CAL.dark1d %ld \n",CAL.dark1d.npixels);
  printf(" npix CAL.errdark1d %ld \n",CAL.errdark1d.npixels);
  printf(" npix CAL.flat1d %ld \n",CAL.flat1d.npixels);
  printf(" npix CAL.errflat1d %ld \n",CAL.errflat1d.npixels);

}


void InitOuputImages(void) {

  SCI.errraw.array=vector_f(SCI.raw.npixels);
  SCI.reduced.array=vector_f(SCI.raw.npixels);
  SCI.errreduced.array=vector_f(SCI.raw.npixels);
  SCI.trimreduced.array=vector_f(SCI.raw.npixels);
  SCI.errtrimreduced.array=vector_f(SCI.raw.npixels);

  SCI.errraw.naxes[0]=SCI.raw.nx;  SCI.errraw.naxes[1]=SCI.raw.ny;                     SCI.errraw.nx=SCI.raw.nx;  SCI.errraw.ny=SCI.raw.ny;                     SCI.errraw.npixels=SCI.raw.npixels;
  SCI.reduced.naxes[0]=SCI.raw.nx;  SCI.reduced.naxes[1]=SCI.raw.ny;                   SCI.reduced.nx=SCI.raw.nx;  SCI.reduced.ny=SCI.raw.ny;                   SCI.reduced.npixels=SCI.raw.npixels;
  SCI.errreduced.naxes[0]=SCI.raw.nx;  SCI.errreduced.naxes[1]=SCI.raw.ny;             SCI.errreduced.nx=SCI.raw.nx;  SCI.errreduced.ny=SCI.raw.ny;             SCI.errreduced.npixels=SCI.raw.npixels;
  SCI.trimreduced.naxes[0]=SCI.raw.nx;  SCI.trimreduced.naxes[1]=SCI.raw.ny;           SCI.trimreduced.nx=SCI.raw.nx;  SCI.trimreduced.ny=SCI.raw.ny;           SCI.trimreduced.npixels=SCI.raw.npixels;
  SCI.errtrimreduced.naxes[0]=SCI.raw.nx;  SCI.errtrimreduced.naxes[1]=SCI.raw.ny;     SCI.errtrimreduced.nx=SCI.raw.nx;  SCI.errtrimreduced.ny=SCI.raw.ny;     SCI.errtrimreduced.npixels=SCI.raw.npixels;

}


void PerformOperations(void) {

  int i,j;
  int nx,ny;
  nx=SCI.raw.nx;
  ny=SCI.raw.ny;

  printf(" Computing reduced and error images\n");

  for(i=0;i<nx;i++) {
    for(j=0;j<ny;j++) {
/*       printf(" i %d j %d RE %f\n",i,j,SCI.reduced.array[i+j*nx]); */
/*       printf("           RA %f\n",SCI.raw.array[i+j*nx]);  */
/*       printf("           DA %f\n",CAL.dark1d.array[i]); */
/*       printf("           FA %f\n",CAL.flat1d.array[i]); */
/*       printf("           ED %f\n",CAL.errdark1d.array[i]); */
/*       printf("           EF %f\n",CAL.errflat1d.array[i]); */
/*       printf("           ER %f\n",SCI.errreduced.array[i+j*nx]); */
      SCI.reduced.array[i+j*nx]=(SCI.raw.array[i+j*nx]-CAL.dark1d.array[i])/CAL.flat1d.array[i];
/*       printf(" i %d j %d RE %f\n",i,j,SCI.reduced.array[i+j*nx]);  */
/*       printf(" Aqui no \n"); */
      SCI.errreduced.array[i+j*nx]=sqrt(((SCI.raw.array[i+j*nx]/CAL.gain + CAL.ron*CAL.ron)+CAL.errdark1d.array[i]*CAL.errdark1d.array[i] + SCI.reduced.array[i+j*nx]*SCI.reduced.array[i+j*nx]*CAL.errflat1d.array[i]*CAL.errflat1d.array[i])/CAL.flat1d.array[i]);
/*       printf(" Aqui si \n"); */
    }
  }
}



void CheckDim(void) {


  if(CAL.dark1d.ny==1 && CAL.errdark1d.ny==1 && CAL.flat1d.ny==1 && CAL.errflat1d.ny==1 );
  else {
    printf(" ERROR: One or more of calibration images are not 1-d FITS files\n");
    exit(1);
  }

  if(CAL.dark1d.nx==CAL.errdark1d.nx && CAL.flat1d.nx==CAL.errflat1d.nx && CAL.dark1d.nx==CAL.flat1d.nx && CAL.dark1d.nx==SCI.raw.nx);
  else {
    printf(" ERROR: X dimension of calibration images do not match with raw image \n");
    exit(1);
  }


}


void SaveOutputImages(void) {

  int status=0;

  printf(" Saving output image %s\n",SCI.trimreduced.file);

  ffinit(&(SCI.trimreduced.filefits),SCI.trimreduced.file,&status);
  if(status) {
    printf(" Couldn't save %s\n",SCI.trimreduced.file);
    fits_report_error(stderr,status);  
    exit(1);
  }
  fits_create_img(SCI.trimreduced.filefits,USHORT_IMG,2,SCI.trimreduced.naxes,&status);
  fits_write_img(SCI.trimreduced.filefits,TFLOAT,1,SCI.trimreduced.naxes[0]*SCI.trimreduced.naxes[1],SCI.trimreduced.array,&status);
  status=0;
  fits_close_file(SCI.trimreduced.filefits,&status);

  printf(" Saving output image %s \n",SCI.errtrimreduced.file);

  ffinit(&(SCI.errtrimreduced.filefits),SCI.errtrimreduced.file,&status);
  if(status) {
    printf(" Couldn't save %s\n",SCI.errtrimreduced.file);
    fits_report_error(stderr,status);  
    exit(1);
  }
  fits_create_img(SCI.errtrimreduced.filefits,USHORT_IMG,2,SCI.errtrimreduced.naxes,&status);
  fits_write_img(SCI.errtrimreduced.filefits,TFLOAT,1,SCI.errtrimreduced.naxes[0]*SCI.errtrimreduced.naxes[1],SCI.errtrimreduced.array,&status);
  status=0;
  fits_close_file(SCI.errtrimreduced.filefits,&status);

  if(status) {
    printf(" There was an error\n");
    exit(1);
  }

}

void TrimImages(void) {

  int ixini,ixfin;
  int iyini,iyfin;
  char st1[200],st2[200],st3[200],st4[200];
  int i,j;
  char *stm,*sto;

  stm=CAL.trimsec;

  
  sto=strstr(stm,"[")+1;
  strcpy(st1,sto);
  strtok(st1,":");
  ixini=atof(st1);
  stm=strstr(sto,":")+1;
  strcpy(st2,stm);
  strtok(st2,",");
  ixfin=atof(st2);
  sto=strstr(stm,",")+1;
  strcpy(st3,sto);
  strtok(st3,":");
  iyini=atof(st3);
  stm=strstr(sto,":")+1;
  strcpy(st4,stm);
  strtok(st4,"]");
  iyfin=atof(st4);

  printf(" Triming [ %d : %d , %d : %d ] images \n",ixini,ixfin,iyini,iyfin);

  SCI.trimreduced.nx=ixfin-ixini+1;
  SCI.trimreduced.ny=iyfin-iyini+1;
  SCI.trimreduced.naxes[0]=ixfin-ixini+1;
  SCI.trimreduced.naxes[1]=iyfin-iyini+1;
  SCI.trimreduced.npixels=SCI.trimreduced.nx*SCI.trimreduced.ny;;
  SCI.errtrimreduced.nx=ixfin-ixini+1;
  SCI.errtrimreduced.ny=iyfin-iyini+1;
  SCI.errtrimreduced.naxes[0]=ixfin-ixini+1;
  SCI.errtrimreduced.naxes[1]=iyfin-iyini+1;
  SCI.errtrimreduced.npixels=SCI.errtrimreduced.nx*SCI.errtrimreduced.ny;;
  
  for(i=ixini-1;i<ixfin-1;i++) {
    for(j=iyini-1;j<iyfin-1;j++) {
      SCI.trimreduced.array[i-ixini+1 + (j-iyini+1)*SCI.trimreduced.nx]=SCI.reduced.array[i+ (j) *SCI.reduced.nx];
      SCI.errtrimreduced.array[i-ixini+1 + (j-iyini+1)*SCI.errtrimreduced.nx]=SCI.errreduced.array[i+ (j) *SCI.errreduced.nx];
    }
  }
  


}

void CloseFiles(void) {

  int status=0;

  fits_close_file(CAL.dark1d.filefits,&status);
  fits_close_file(CAL.errdark1d.filefits,&status);
  fits_close_file(CAL.flat1d.filefits,&status);
  fits_close_file(CAL.errflat1d.filefits,&status);
  fits_close_file(SCI.raw.filefits,&status);
  fits_close_file(SCI.errraw.filefits,&status);

}

