#include "modulos.h"


int main()
{
  
  char inputfile[100]="";
  int status=0;
  int nfound, anynull;
  fitsfile *inputimage;
  long naxes[2], fpixel,  npixels;
  float datamin, datamax, nullval;
  
  float *array;

  int bitpix;
  int i;


  char outputfile[100]="";
  fitsfile *outputimage;
  float tr[6];

  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;

  
  /* Leo la imagen FITS de entrada */
  
  printf(" Input image: ");
  reads(inputfile,inputfile);

  printf(" BITPIX for output image (8,16,32,-32,64) : ");
  bitpix=readi(-32);

  printf(" Output image: ");
  reads(outputfile,outputfile);


  if( ffopen(&inputimage, inputfile, READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(inputimage, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  npixels=naxes[0]*naxes[1];
  if((array=malloc(npixels*sizeof(float)))==NULL) printf("I cannot dimension buffer of %ld elements\n",npixels);
  fpixel=1;
  nullval=0;
  datamin=1.0e30;
  datamax=-1.0e30;
  printf("...Reading image %s \n",inputfile);
  if(fits_read_img(inputimage, TFLOAT, fpixel, npixels, &nullval, array, &anynull, &status )) fits_report_error(stderr,status);


  for(i=0;i<npixels;i++) printf("%d  %f \n",i,array[i]); 

  cpgopen("?");
  cpgwnad(0.,naxes[0],0.,naxes[1]);
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  cpggray(array,naxes[0],naxes[1],1,naxes[0],1,naxes[1],40,-10,tr);


  ffinit(&outputimage,outputfile,&status);
/*   fits_delete_key(outputimage,"BITPIX",&status); */
  fits_create_img(outputimage,bitpix,2,naxes,&status);
/*   ffcphd(inputimage,outputimage,&status); */
  fits_write_img(outputimage,TFLOAT,1,naxes[0]*naxes[1],array,&status);


  fits_close_file(outputimage,&status);
  fits_close_file(inputimage,&status);
  
  if(status) fits_report_error(stderr,status);


  return 0;
}

