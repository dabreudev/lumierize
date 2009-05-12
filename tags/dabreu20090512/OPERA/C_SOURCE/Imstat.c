#include "modulos.h"


int main()
{

  float mean,sigma;
  float first,median,third;
  float xbmin,xbmax,ybmin,ybmax;
  float xfmin,xfmax,yfmin,yfmax;
  char opt[10];
  char cnul;
  char tittle[200];
  float imin,imax,jmin,jmax;
  /* Varaibales para leer el FITS */
  char inputfile[251];
  int status=0;
  int nfound, anynull;
  fitsfile *image;
  long naxes[2], fpixel, npixels, ii,jj;
  float datamin, datamax, nullval;
  float *buffer;
  /* Variables del dibujo */
/*   float factorarr=1; */
  float tr[6];
  /* Variables para la fotometria */
  float sky, flux,object;
  int npixobj, npixsky;
  int iobj,nobj;
  float mag[1000],cts[1000],logcts[1000],a,b;
  float col[1000];
  int log[1000];
  int i;
  float ctsmin,ctsmax,magmin,magmax; 
  FILE *statfile;
  char filestat[51];
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;

  fflush(NULL);
  cpgbeg(0,"?",1,1);

  printf(" Input FITS file: ");
  reads("",inputfile);
  
    
/* Leo la imagen FITS de entrada */
  if( ffopen(&image, inputfile, READONLY, &status)) fits_report_error(stderr,status);
  if(fits_read_keys_lng(image, "NAXIS", 1, 2, naxes, &nfound, &status)) fits_report_error(stderr,status);
  npixels=naxes[0]*naxes[1];
  if((buffer=malloc(npixels*sizeof(float)))==NULL) printf("I cannot dimension buffer of %ld elements\n",npixels);
  fpixel=1;
  nullval=0;
  datamin=1.0e30;
  datamax=-1.0e30;
  printf("...Reading image %s \n",inputfile);
  if(fits_read_img(image, TFLOAT, fpixel, npixels, &nullval, buffer, &anynull, &status )) fits_report_error(stderr,status);




  printf(" Name of output file with statistics: ");
  scanf("%s",filestat);
  statfile=fopen(filestat,"w");
  fprintf(statfile,"#Image    Npix Mean  Stddev  Min Max Mode  1º Quartil 3º Quartil\n"); 

  printf("...Computing statistics \n");
  MinMax(npixels,buffer,&datamin,&datamax);
  

  mean=StMedia(npixels,buffer,&sigma);
  Quartil(npixels,buffer,&first,&median,&third);
  fprintf(statfile," %s     %ld  %f   %f    %f   %f   %f   %f   %f\n",inputfile,npixels,mean,sigma,datamin,datamax,median,first,third);

  /* Dibujo la imagen  */


  cpgscf(2);
  cpgask(0);
  //cpglab("pixel","pixel","");
  cpgvstd();
  //  cpgswin(0.,naxes[0],0.,naxes[1]);

  cpgwnad(0.,naxes[0],0.,naxes[1]);

  
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  cpggray(buffer,naxes[0],naxes[1],1,naxes[0],1,naxes[1],datamax,datamin,tr);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  
  xbmin=1;xbmax=naxes[0];ybmin=1;ybmax=naxes[1];

 

    /* Menu principal */
  exit(0);

  while(opt[0]!='E') {
    
    printf("\n\n  Z Zoom image with cursor\n");
    printf("  O Zoom out\n");
    printf("  F Fotometry (number of counts)\n");
    printf("  U Calibration curve\n");
    printf("  R Read photometry file\n");
    printf("  W Write photometry file\n");
    printf("  C Change cuts max and min\n");
    printf("  E Exit\n");
    fflush(NULL);
    fflush(stdin);
    fflush(stdout);
    scanf("%s",opt);
    fflush(NULL);
    fflush(stdin);
    fflush(stdout);


    switch (opt[0]) {
    case 'Z' :
      cpgsci(2);
      printf(" Press bottom left square with mouse...\n");
      cpgcurs(&xbmin,&ybmin,&cnul);
      printf(" Press top right square with mouse...\n");
      cpgband(2,1,xbmin,ybmin,&xbmax,&ybmax,&cnul);
      cpgsci(1);
      break;
    case 'O' : 
      xbmin=1;xbmax=naxes[0];ybmin=1;ybmax=naxes[1];
      break;

    case 'C' :

      
      printf(" Current Maximum: %f Minimum %f \n",datamax,datamin);
      printf(" New minimum: ");
      scanf("%f",&datamin);
      printf(" New maximum: ");
      scanf("%f",&datamax);
      break;
      
    case 'F' :
      cpgsci(2);
      printf(" ####################\n Region over the object\n");
      printf(" Press bottom left square of rectangle to compute photometry...\n");
      cpgcurs(&xfmin,&yfmin,&cnul);
      printf(" Press top right square of rectangle to compute photometry...\n");
      cpgband(2,1,xfmin,yfmin,&xfmax,&yfmax,&cnul);
      cpgsci(1);
      imin=(int)xfmin;imax=(int)xfmax;jmin=(int)yfmin;jmax=(int)yfmax;
      
      flux=0;
      npixobj=0;
      for (ii=imin ;ii<imax+1;ii++) {
	for (jj=jmin ;jj<jmax+1;jj++) {
	  flux+=buffer[ii+naxes[0]*jj];
	  npixobj++;
	}
      }


      printf(" Object rectangle  [%d:%d,%d:%d]\n",(int)imin,(int)imax,(int)jmin,(int)jmax);

      printf(" ####################\n Region over the sky\n");
      printf(" Press bottom left square of rectangle to compute sky...\n");
      cpgcurs(&xfmin,&yfmin,&cnul);
      printf(" Press top right square of rectangle to compute sky...\n");
      cpgband(2,1,xfmin,yfmin,&xfmax,&yfmax,&cnul);
      cpgsci(1);
      imin=(int)xfmin;imax=(int)xfmax;jmin=(int)yfmin;jmax=(int)yfmax;
      sky=0;
      npixsky=0;
      for (ii=imin ;ii<imax+1;ii++) {
	for (jj=jmin ;jj<jmax+1;jj++) {
	  sky+=buffer[ii+naxes[0]*jj];
	  npixsky++;
	}
      }

      printf(" Sky rectangle  [%d:%d,%d:%d]\n",(int)imin,(int)imax,(int)jmin,(int)jmax);

      object=flux-npixobj*sky/npixsky;
      MCP1(nobj,logcts,mag,&a,&b);
      printf(" ####################\n");
      
      printf(" Total flux (object+sky): %f\n",flux);
      printf(" Object flux: %f\n",object);
      printf(" Sky flux / pixel: %f\n",sky/npixsky);
      printf(" Number of pixels from the object: %d\n",npixobj);
      printf(" Number of pixels from the sky: %d\n",npixsky);
      //printf("Photometric calibration mag= %f + %f log10(counts)",a,b); 
      printf(" Magnitude with current calibration: %f\n",a+b*log10(object));
      printf(" ####################\n");
      //      printf(" Mean: %f\n",flux/npixobj);
      printf(" Input magnitude of  object %d: ",nobj+1);
      scanf("%f",mag+nobj);

      cts[nobj]=object;
      logcts[nobj]=log10(object);
      nobj++;

      break;
    case 'E' :
      exit(0);
      break;
    case 'U' :
      printf(" Calibration curve for the image using %d objects\n",nobj);
      cpgpage();
      pgLimits(nobj,logcts,&ctsmin,&ctsmax);
      pgLimits(nobj,mag,&magmin,&magmax);
      
      //printf(" CTS min %f Max %f \n",ctsmin,ctsmax);
      //printf(" MAG min %f Max %f \n",magmin,magmax);
      cpgswin(ctsmin,ctsmax,magmin,magmax);
      cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);  
      cpgsch(1.);
      MCP1(nobj,logcts,mag,&a,&b);
      sprintf(tittle,"Photometric calibration mag= %f + %f log10(counts)",a,b); 
      cpglab("log (object counts)","Magnitude",tittle);
      cpgsch(3.);
     for(iobj=0;iobj<nobj;iobj++) {
	printf(" Object %d   magnitude %f  counts %f  log(counts) %f \n",iobj+1,mag[iobj],cts[iobj],logcts[iobj]);
	cpgpt1(logcts[iobj],mag[iobj],3);
      }

      cpgsch(1.);
      printf(" Calibration m= %f + %f * log10(object counts)\n",a,b);
      cpgmove(ctsmin,a+b*ctsmin);
      cpgdraw(ctsmax,a+b*ctsmax);
      cpgcurs(&nullval,&nullval,&cnul);
      break;
    case 'W':
      printf(" Name of file with fotometry: ");
      scanf("%s",filestat);
      statfile=fopen(filestat,"w");
      fprintf(statfile,"#Object    Counts    Mag\n"); 
      for(iobj=0;iobj<nobj;iobj++) 
	fprintf(statfile," %d   %f  %f\n",iobj+1,cts[iobj],mag[iobj]);
      fclose(statfile);
      break;
    case 'R' :
      for(i=0;i<100;i++) log[i]=0;      
      printf(" Name of file with fotometry: ");
      scanf("%s",filestat);
      ReadNumcol(filestat,2,col,log,&nobj);
      iobj=0;
      for(i=0;i<nobj;i++) {
	if(log[i]) {
	  cts[iobj]=col[i];
	  logcts[iobj]=log10(col[i]);
	  iobj++; 

	}


      }
      iobj=0;
      ReadNumcol(filestat,3,col,log,&nobj);
      for(i=0;i<nobj;i++) {
	if(log[i]) {
	  mag[iobj]=col[i];
	  printf("  N.  %d  Cts %f  Mag %f\n",iobj+1,cts[iobj],mag[iobj]); 
	  iobj++;
	}
	


      }
      nobj=iobj;
   


    }

    /*

    datamin=1.0e30;
    datamax=-1.0e30;
    for (ii=xbmin ;ii<xbmax+1;ii++) {
      for (jj=ybmin ;jj<ybmax+1;jj++) {
	
	if( buffer[ii+naxes[0]*jj]< datamin) datamin=buffer[ii+naxes[0]*jj];
	if(buffer[ii+naxes[0]*jj]> datamax) datamax=buffer[ii+naxes[0]*jj];
      }
    }
    */
    cpgpage();
    cpgvstd();
    cpgwnad(xbmin,xbmax,ybmin,ybmax);
    //cpgswin(xbmin,xbmax,ybmin,ybmax);

    printf(" Drawing subimage [%f:%f,%f:%f]\n",xbmin,xbmax,ybmin,ybmax);

    cpggray(buffer,naxes[0],naxes[1],xbmin,xbmax,ybmin,ybmax,datamax,datamin,tr);
    cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);    
  }

  return 0;
}
