#include "modulos.h"

void ApertureFot(float *ima, int nx, int ny,float xc, float yc,float rkron,float sky);

int main()
{

  float mean,sigma;
  float first,median,third;
  float xbmin,xbmax,ybmin,ybmax;
  float xfmin,xfmax,yfmin,yfmax;
  float xcenter,ycenter,radius,rkron;
  float xradius,yradius;
  char opt='E';
  char cnul;
  char tittle[200];
  int imin,imax,jmin,jmax;
  /* Varaibales para leer el FITS */
  char inputfile[51];
  int status=0;
  int nfound, anynull;
  fitsfile *image;
  long naxes[2], fpixel, npixels, ii,jj;
  float datamin, datamax, nullval;
  float *buffer;
  float *bufferstat;
  /* Variables del dibujo */
/*   float factorarr=1; */
  float tr[6];
  /* Variables para la fotometria */
  float sky=0, flux,object;
/*   float  noise, snr; */
  int npixobj, npixsky;
  int iobj,nobj;
  float mag[1000],cts[1000],logcts[1000],a,b;
  float col[1000];
  int log[1000];
  int i;
  float ctsmin,ctsmax,magmin,magmax; 
  FILE *fotfile;
  char filefot[51]="/dev/null";
  char statfile[51]="/dev/null";
  FILE *fst;
  int plotflag=0;
  
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
  printf("...Computing minimum and maximum cuts \n");
/*   //  for (ii=0 ;ii<npixels;ii++) { */
/*     //printf(" buffer %d %e \n",ii,buffer[ii]); */
/*   //  if( buffer[ii]< datamin) datamin=buffer[ii]; */
/*   //  if(buffer[ii]> datamax) datamax=buffer[ii]; */
/*   //} */
  
  mean=StMedia(npixels,buffer,&sigma);
  datamin=mean-sigma*3;
  datamax=mean+sigma*3;
  printf(" Datamin cut %f Datamax cut %f \n",datamin,datamax);


  /* Dibujo la imagen  */


  cpgscf(2);
  cpgask(0);
/*   //cpglab("pixel","pixel",""); */
/*   //cpgswin(0.,naxes[0],0.,naxes[1]); */
/*   //cpgsvp(0.1,0.9,0.1,0.9); */
  cpgvstd();
  cpgwnad(0.,naxes[0],0.,naxes[1]);

  
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  cpggray(buffer,naxes[0],naxes[1],1,naxes[0],1,naxes[1],datamax,datamin,tr);
  cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);
  
  xbmin=1;xbmax=naxes[0];ybmin=1;ybmax=naxes[1];

 
  nobj=0;
    /* Menu principal */

  do {
    
    printf("\n\n  Z Zoom image with cursor\n");
    printf("  O Zoom out\n");
    printf("  S Print rough statistics for this image\n");
    printf("  P Print rough statistics for a given region\n");
    printf("  F Photometry (number of counts) in reactangular area\n");
    printf("  C Photometry (number of counts) in circular aperture\n");
    printf("  G Grow curve in aperture photometry\n");    
    printf("  U Calibration curve\n");
    printf("  R Read photometry file\n");
    printf("  W Write photometry file\n");
    printf("  M Change cuts max and min\n");
    printf("  D Draw image in 3D\n");
    printf("  B Draw image in 3D with bars\n");
    printf("  I Draw image in a gray map (default)\n");
    printf("  E Exit\n");
    fflush(NULL);
    fflush(stdin);
    fflush(stdout);
    opt=readc(opt);
    fflush(NULL);
    fflush(stdin);
    fflush(stdout);


    switch (opt) {
    case 'S':
    case 's':
      printf(" Name of output file with statistics: ");
      reads(statfile,statfile);
      mean=StMedia(npixels,buffer,&sigma);
      Quartil(npixels,buffer,&first,&median,&third);
      if((fst=fopen(statfile,"w"))==NULL) printf(" An error ocurred while opening %s\n",statfile);
      fprintf(fst,"#     Image                               mean    sigma   first_quartil   median   third_quartil\n");
      fprintf(fst,"%-51s   %10.7f  %10.7f  %10.7f  %10.7f  %10.7f\n",inputfile,mean,sigma,first,median,third);
      printf("#     Image                               mean    sigma   first_quartil   median   third_quartil\n"); 
      printf("%-51s   %10.7f  %10.7f  %10.7f  %10.7f  %10.7f\n",inputfile,mean,sigma,first,median,third);
      fclose(fst);
      break;
    case 'Z' :
    case 'z' :
      cpgsci(2);
      printf(" Press bottom left square with mouse...\n");
      cpgcurs(&xbmin,&ybmin,&cnul);
      printf(" Press top right square with mouse...\n");
      cpgband(2,1,xbmin,ybmin,&xbmax,&ybmax,&cnul);
      cpgsci(1);
      break;
    case 'O' : 
    case '0' : 
      xbmin=1;xbmax=naxes[0];ybmin=1;ybmax=naxes[1];
      break;

    case 'M' :
    case 'm' :

      
      printf(" Current Maximum: %f Minimum %f \n",datamax,datamin);
      printf(" New minimum: ");
      datamin=readf(datamin);
      printf(" New maximum: ");
      datamax=readf(datamax);
      break;
      
    case 'F' :
    case 'f' :
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
/*       //printf("Photometric calibration mag= %f + %f log10(counts)",a,b);  */
      printf(" Magnitude with current calibration: %f\n",a+b*log10(object));
      printf(" ####################\n");
/*       //      printf(" Mean: %f\n",flux/npixobj); */
      printf(" Input magnitude of  object %d: ",nobj+1);
      mag[nobj]=readf(mag[nobj]);

      cts[nobj]=object;
      logcts[nobj]=log10(object);
      nobj++;

      break;
    case 'P' :
    case 'p' :
      cpgsci(2);
      printf(" Press bottom left square of rectangle to compute statistics...\n");
      cpgcurs(&xfmin,&yfmin,&cnul);
      printf(" Press top right square of rectangle to compute statistics...\n");
      cpgband(2,1,xfmin,yfmin,&xfmax,&yfmax,&cnul);
      cpgsci(1);
      imin=(int)xfmin;imax=(int)xfmax;jmin=(int)yfmin;jmax=(int)yfmax;
      sky=0;
      npixsky=0;
      bufferstat=vector_f((imax-imin+1)*(jmax-jmin+1));
      printf(" Alloc\n");
      for (ii=imin ;ii<imax+1;ii++) {
	for (jj=jmin ;jj<jmax+1;jj++) {
	  bufferstat[npixsky]=buffer[ii+naxes[0]*jj];
	  npixsky++;
	}
      }
      mean=StMedia(npixsky,bufferstat,&sigma);
      Quartil(npixsky,bufferstat,&first,&median,&third); 
      printf("#     Image                               mean    sigma   first_quartil   median   third_quartil\n"); 
      printf("%-51s[%d:%d,%d:%d]   %10.7f  %10.7f  %10.7f  %10.7f  %10.7f\n",inputfile,(int)imin,(int)imax,(int)jmin,(int)jmax,mean,sigma,first,median,third);
      printf(" Sigma with quartil: %f\n",(third-first)/1.35);

      free(bufferstat);
      break;
    case 'C':
    case 'c':
      cpgsci(2);
      printf(" ####################\n Region over the object\n");
      printf(" Press center of the aperture...\n");
      cpgcurs(&xcenter,&ycenter,&cnul);
/*       //printf(" Object center ads sa : %f , %f  ",xcenter,ycenter); */

      printf(" Press up to the radius to compute photometry...\n");
      cpgband(1,1,xcenter,ycenter,&xradius,&yradius,&cnul);
      radius=sqrt((xradius-xcenter)*(xradius-xcenter)+(yradius-ycenter)*(yradius-ycenter));
/*       //printf("  Radius %f\n",radius); */
      cpgsci(1);
/*       //printf(" Antes de entrar nx %d ny %d\n",(int)naxes[0],(int)naxes[1]); */
      flux=circ_aper(buffer,naxes[0],naxes[1],xcenter,ycenter,radius,&npixobj);
      printf(" Con la circular %f\n",flux);
/*       flux=elip_aper(buffer,(int)naxes[0],(int)naxes[1],xcenter,ycenter,radius,radius,0,&npixobj); */
/*       printf(" Con la eliptica %f\n",flux); */
/*       flux=elip_aper(buffer,(int)naxes[0],(int)naxes[1],xcenter,ycenter,radius/2.,radius*2.,3.14,&npixobj); */
/*       printf(" Con la eliptica 2 %f\n",flux); */
/*       flux=elip_aper(buffer,(int)naxes[0],(int)naxes[1],xcenter,ycenter,radius*2.,radius/2.,3.14,&npixobj); */
/*       printf(" Con la eliptica 3 %f\n",flux); */

      printf(" Object center  %d , %d  Radius %f\n",(int)xcenter,(int)ycenter,radius);
      
      printf(" ####################\n Region over the sky\n");
      printf(" Press center of the aperture to compute sky...\n");
      cpgcurs(&xcenter,&ycenter,&cnul);
      printf(" Press up to the radius to compute photometry...\n");
      cpgband(1,1,xcenter,ycenter,&xradius,&yradius,&cnul);
      radius=sqrt((xradius-xcenter)*(xradius-xcenter)+(yradius-ycenter)*(yradius-ycenter));
      cpgsci(1);
      sky=circ_aper(buffer,naxes[0],naxes[1],xcenter,ycenter,radius,&npixsky);
      
      printf(" Sky center  %d , %d  Radius %f\n",(int)xcenter,(int)ycenter,radius);
      
      object=flux-npixobj*sky/npixsky;
      MCP1(nobj,logcts,mag,&a,&b);
      printf(" ####################\n");
      
      printf(" Total flux (object+sky): %f\n",flux);
      printf(" Object flux: %f\n",object);
      printf(" Sky flux / pixel: %f\n",sky/npixsky);
      printf(" Number of pixels from the object: %d\n",npixobj);
      printf(" Number of pixels from the sky: %d\n",npixsky);
/*       //printf("Photometric calibration mag= %f + %f log10(counts)",a,b);  */
      printf(" Magnitude with current calibration: %f\n",a+b*log10(object));
      printf(" ####################\n");
/*       //      printf(" Mean: %f\n",flux/npixobj); */
      printf(" Input magnitude of  object %d: ",nobj+1);
      mag[nobj]=readf(mag[nobj]);
      
      cts[nobj]=object;
      logcts[nobj]=log10(object);
      nobj++;
      break;
    case 'G':
    case 'g':
      printf(" ####################\n Press center of the aperture...\n");
      cpgcurs(&xcenter,&ycenter,&cnul);
      rkron=sqrt(xcenter*xcenter+ycenter*ycenter);
      ApertureFot(buffer,naxes[0],naxes[1],xcenter,ycenter,rkron,sky);

      break;
    case 'E' :
    case 'e' :
/*       exit(0); */
      break;
    case 'U' :
    case 'u' :
      printf(" Calibration curve for the image using %d objects\n",nobj);
      cpgpage();
      pgLimits(nobj,logcts,&ctsmin,&ctsmax);
      pgLimits(nobj,mag,&magmin,&magmax);
      
/*       //printf(" CTS min %f Max %f \n",ctsmin,ctsmax); */
/*       //printf(" MAG min %f Max %f \n",magmin,magmax); */
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
    case 'w':
      printf(" Name of file with fotometry: ");
      reads(filefot,filefot);
      fotfile=fopen(filefot,"w");
      fprintf(fotfile,"#Object    Counts    Mag\n"); 
      for(iobj=0;iobj<nobj;iobj++) 
	fprintf(fotfile," %d   %f  %f\n",iobj+1,cts[iobj],mag[iobj]);
      fclose(fotfile);
      break;
    case 'D':
    case 'd':
      plotflag=2;
      break;
    case 'B':
    case 'b':
      plotflag=1;
      break;
    case 'I':
    case 'i':
      plotflag=0;
      break;
    case 'R' :
    case 'r' :
      for(i=0;i<100;i++) log[i]=0;      
      printf(" Name of file with fotometry: ");
      reads(filefot,filefot);
      ReadNumcol(filefot,2,col,log,&nobj);
      iobj=0;
      for(i=0;i<nobj;i++) {
	if(log[i]) {
	  cts[iobj]=col[i];
	  logcts[iobj]=log10(col[i]);
	  iobj++; 

	}


      }
      iobj=0;
      ReadNumcol(filefot,3,col,log,&nobj);
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
/*     //    cpgwnad(xbmin,xbmax,ybmin,ybmax); */
    

    if(plotflag==2) {
      cpgvstd();
      cpgwnad(0,naxes[0],0,naxes[1]);
      cpgdraw3d(buffer,naxes[0],naxes[1],(float)naxes[0],30.);
    }
    else  if(plotflag==1) {
      cpghist3d(buffer,naxes[0],naxes[1],xbmin,xbmax,ybmin,ybmax);
    }
    else {
      cpgvstd();
      cpgwnad(xbmin,xbmax,ybmin,ybmax);
      printf(" Drawing subimage [%f:%f,%f:%f]\n",xbmin,xbmax,ybmin,ybmax);
      cpggray(buffer,naxes[0],naxes[1],xbmin,xbmax,ybmin,ybmax,datamax,datamin,tr);
      cpgbox("BCTNS",0.0,0,"BCTNS",0.0,0);    
    }
  }  while(opt!='E' && opt!='e');
  cpgclos();  
  cpgend();    
  return(0); 
  
}









void ApertureFot(float *ima, int nx, int ny,float xc, float yc,float rkron,float sky)
{

/*   //Subrutina para calcular fotometria de apertura */
/*   // De momento lo hago copiando lo de ventana en otra matriz de */
/*   // modo que se necesita el doble de memoria. */
  float tr[6]={0.,1.,0.,0.,0.,1.};
  int pg;
  float r;
  int npix;
  int i;
  float fnul;
  char cnul;
  float apfot;
  cpgqid(&pg);
  cpgopen("/xserve");
  cpgsvp(0.1,0.4,0.1,0.9);
  cpgwnad((int)xc-30,(int)xc+30,1.,(int)yc+50);
  cpggray(ima,nx,ny,(int)xc-30,(int)xc+30,1,(int)yc+50,sky-1.5*sqrt(sky),sky+3.5*sqrt(sky),tr);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  cpgsci(2);
  cpgslw(3);
  cpgsls(2);
  cpgslw(1);
  cpgsls(1);
  cpgsci(1);
  cpgsvp(0.5,0.9,0.1,0.9);
  cpgswin(0.,rkron*5,0.,sky*25);
/*   //printf(" Aqui no se quen\n"); */
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
/*   //  printf(" Seque da\n"); */
  cpgsch(2.);
/*   //printf("Entra\n"); */
 
  cpgsfs(2);
  for(i=0;i< 15;i++) {
    cpgsci(i);
    r=rkron*5/15*i;
/*     //printf(" IMAAAA 590 21 %f\n",ima[590+nx*21]); */
    apfot=circ_aper(ima,nx,ny,xc,yc,r,&npix);
    cpgsvp(0.1,0.4,0.1,0.9);
    cpgwnad((int)xc-30,(int)xc+30,1.,(int)yc+50);
    cpgcirc(xc,yc,r);
    cpgsvp(0.5,0.9,0.1,0.9);
    cpgswin(0.,rkron*5,0.,sky*25);
    printf(" radio %f fot %f fot-sky %f npix %d\n",r,apfot,apfot-npix*sky,npix);
    cpgpt1(r,apfot-npix*sky,3);
  }

  printf(" Los pixeles cnetrle s %f %f SKY %f\n",xc,yc,sky);
  printf(" Rkron %f\n",rkron);
  cpgcurs(&fnul,&fnul,&cnul);
  cpgclos();
  cpgslct(pg);
}








