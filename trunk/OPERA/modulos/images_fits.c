#include "modulos.h" 


void ReadImage(struct image *im) {
 
  int status=0; 
  int nfound, anynull;
  long fpixel;
  int nullval;

  if(ffopen(&((*im).filefits), (*im).file, READONLY, &status)) fits_report_error(stderr,status);  
  if(fits_read_keys_lng((*im).filefits, "NAXIS", 1, 2, (*im).naxes, &nfound, &status)) fits_report_error(stderr,status);

  if(nfound==1)  {
    (*im).npixels=(*im).naxes[0];
    (*im).naxes[1]=1;
  }
  else (*im).npixels=(*im).naxes[0]*(*im).naxes[1];
  if((*im).aloc_flag) free((*im).array);
  (*im).array=vector_f((*im).npixels);
  (*im).aloc_flag=1;
  fpixel=1;
  nullval=0;
  (*im).datamin=0;
  (*im).datamax=0;
  (*im).nx=(*im).naxes[0];
  (*im).ny=(*im).naxes[1];
  printf("...Reading image %s \n",(*im).file);
  if(fits_read_img((*im).filefits, TFLOAT, fpixel, (*im).npixels, &nullval, (*im).array, &anynull, &status )) fits_report_error(stderr,status);
  
  ffclos((*im).filefits,&status);
  if(status) {
    fits_report_error(stderr,status);
    exit(status);
  }

}
void ReadWCSImage(struct wcsimage *im) {
 
  int status=0; 
  int nfound, anynull;
  long fpixel;
  int nullval;
  fitsfile *wcstempheader;
  int lhead,nbfits;
  char *header;
  char tempfile[200];
  long rand_i;
  int time_i;

  printf(" Antes ffopen %s\n",(*im).image.file);


  
  if(ffopen(&((*im).image.filefits), (*im).image.file, READONLY, &status)) fits_report_error(stderr,status);  
  printf(" Despue ssasd a\n");
  if(fits_read_keys_lng((*im).image.filefits, "NAXIS", 1, 2, (*im).image.naxes, &nfound, &status)) fits_report_error(stderr,status);

  printf(" Despues ffopen\n");

  if(nfound==1)  {
    (*im).image.npixels=(*im).image.naxes[0];
    (*im).image.naxes[1]=1;
  }
  else (*im).image.npixels=(*im).image.naxes[0]*(*im).image.naxes[1];
  if((*im).image.aloc_flag) free((*im).image.array);
  (*im).image.array=vector_f((*im).image.npixels);
  (*im).image.aloc_flag=1;
  fpixel=1;
  nullval=0;
  (*im).image.datamin=0;
  (*im).image.datamax=0;
  (*im).image.nx=(*im).image.naxes[0];
  (*im).image.ny=(*im).image.naxes[1];
  printf("...Reading image %s \n",(*im).image.file);
  if(fits_read_img((*im).image.filefits, TFLOAT, fpixel, (*im).image.npixels, &nullval, (*im).image.array, &anynull, &status )) fits_report_error(stderr,status);
  rand_i=random();
  time_i=(int)time(NULL);
  sprintf(tempfile,"temp%010d_%012ld_header.fits",time_i,rand_i);
  unlink(tempfile);
  fprintf(stderr," Using %s as temp file \n",tempfile); 
  ffinit(&wcstempheader,tempfile,&status);
  fits_copy_header((*im).image.filefits,wcstempheader,&status);
  ffclos(wcstempheader,&status);
  if ((header = fitsrhead (tempfile, &lhead, &nbfits)) == NULL) {
    printf(" Couldn`t open %s\n",tempfile);
    exit(1);
  }
  if(((*im).wcs=wcsinit(header))==NULL) {
    printf(" Warning: %s does not contain WCS information\n",tempfile);
    (*im).wcsset=0;
  }
  else (*im).wcsset=1;
  unlink(tempfile);
  
  ffclos((*im).image.filefits,&status);
  if(status) {
    fits_report_error(stderr,status);
    exit(status);
  }
  free(header); 
}

void SaveImage(struct image im) {
 
  int status=0; 
  long naxes[2];


  ffinit(&(im.filefits),im.file,&status);
  if(status) fits_report_error(stderr,status);
  naxes[0]=im.nx;
  naxes[1]=im.ny;
  fits_create_img(im.filefits,-32,2,naxes,&status);
  printf("...Saving image %s \n",im.file);
  fits_write_img(im.filefits,TFLOAT,1,im.nx*im.ny,im.array,&status);
  fits_close_file(im.filefits,&status);
  if(status) {
    fits_report_error(stderr,status);
    exit(status);
  }
}

void SaveImage_16(struct image im) {
 
  int status=0; 
  long naxes[2];


  ffinit(&(im.filefits),im.file,&status);
  if(status) fits_report_error(stderr,status);
  naxes[0]=im.nx;
  naxes[1]=im.ny;
  fits_create_img(im.filefits,16,2,naxes,&status);
  printf("...Saving image %s \n",im.file);
  fits_write_img(im.filefits,TFLOAT,1,im.nx*im.ny,im.array,&status);
  fits_close_file(im.filefits,&status);
  if(status) {
    fits_report_error(stderr,status);
    exit(status);
  }

}



void PlotImage(struct image *im) {

  float mean,sigma;
  float first,median,third;
  float tr[6];
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;


  if(im->datamin==0 && im->datamax==0) {
    mean=StMedia(im->npixels,im->array,&sigma); 
    im->datamin=mean-sigma*3;
    im->datamax=mean+sigma*3;
    
    Quartil(im->npixels,im->array,&first,&median,&third);
    im->datamin=median-(median-first)*3;
    im->datamax=median+(third-median)*6;
    
    printf(" mean %f sigma %f\n",mean,sigma);
    printf(" first %f third %f\n",first,third);
  }


  cpgscf(2);
  cpgask(0);
  cpgeras();
  cpgvstd();
  cpgwnad(0.,im->naxes[0],0.,im->naxes[1]);
  cpggray(im->array,im->naxes[0],im->naxes[1],1,im->naxes[0],1,im->naxes[1],im->datamax,im->datamin,tr);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);

}

void PlotImage_sec(struct image *im, int x1, int x2, int y1, int y2) {

  float mean,sigma;
  float first,median,third;
  float tr[6];
  float *buffer;
  int npix;
  int i,j,k;
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;

  npix=(x2-x1+1)*(y2-y1+1);
  buffer=vector_f(npix);

  k=0;
  for(i=x1;i<=x2;i++) {
    for(j=y1;j<=y2;j++) {
      buffer[k]=im->array[i+j*im->nx];
      k++;
    }
  }
  Quartil(npix,buffer,&first,&median,&third);

  cpgscf(2);
  cpgask(0);
  cpgeras();
  cpgvstd();
  cpgwnad(x1,x2,y1,y2);
  cpggray(im->array,im->naxes[0],im->naxes[1],x1,x2,y1,y2,median+(median-first)*8,median-(median-first)*3,tr);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  free(buffer);

}




void PlotImage_sec_cuts(struct image *im, int x1, int x2, int y1, int y2, float bg, float fg) {

  float tr[6];
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;


  cpgscf(2);
  cpgask(0);
  cpgeras();
  cpgvstd();
  cpgwnad(x1,x2,y1,y2);
  cpggray(im->array,im->naxes[0],im->naxes[1],x1,x2,y1,y2,fg,bg,tr);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);

}

