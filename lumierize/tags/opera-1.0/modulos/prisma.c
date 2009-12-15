#include "modulos.h"

#define DEBUG 0


float pr_ldo2pix(float ldo,struct disper_prism DP)
{
  return((DP.A+DP.B/(ldo-DP.C))/DP.tampix);  
/*   return(1.5*ldo+0.5);  */
}

float pr_pix2ldo(float pix,struct disper_prism DP)
{
  return(DP.C-DP.B/(DP.A-pix*DP.tampix));  
/*   return(pix/1.5-.5);  */
}

float pr_dpdl(float ldo,struct disper_prism DP)
{
  return(-DP.B/((ldo-DP.C)*(ldo-DP.C))/DP.tampix);  
/*   return(1.5);  */
}


float pr_dldp(float pix,struct disper_prism DP)
{
  return(-DP.B/((DP.A-pix*DP.tampix)*(DP.A-pix*DP.tampix))*DP.tampix);  
/*   return(1/1.5);  */
}


void SP_pix2ldo(struct spectrum *SP_lam, struct spectrum SP_pix, struct disper_prism DP) {
  
  /* Paso de un espectro por pixel a un espectro por lambda */
  /* De paso alocatea el espectro de lambda */

  int i;
  float pixreal1;
  float pixreal2;
  float ldo1,ldo2;
  float ldomin,ldomax,deltaldo;
/*   float *pixels; */
/*   static int npix_buf=0; */

  if(DEBUG) printf(" Estoy en esta zona de aqui\n");

  /*   if(npix_buf!=SP_pix.nx) { */
  /*     if(npix_buf!=0) free(pixels); */
  /*     pixels=vector_f(SP_pix.nx); */
  /*     for(i=0;i<SP_pix.nx;i++) pixels[i]=(float)(i); */
  /*   } */
  ldomin=pr_pix2ldo(1.,DP);
  ldomax=pr_pix2ldo(SP_pix.nx,DP);
  if(DP.B>0) SWAP(ldomin,ldomax);
  deltaldo=(ldomax-ldomin)/(SP_pix.nx-1.);
  (*SP_lam).ldomin=ldomin;
  (*SP_lam).deltaldo=deltaldo;
  if((*SP_lam).aloc_flag==1) free((*SP_lam).spec);
  (*SP_lam).spec=vector_f(SP_pix.nx);
  (*SP_lam).nx=SP_pix.nx;
  (*SP_lam).naxes[0]=SP_pix.nx;
  (*SP_lam).npixels=SP_pix.nx;
  (*SP_lam).aloc_flag=1;
  if((*SP_lam).alocldo_flag==1) free((*SP_lam).ldo);
  (*SP_lam).ldo=vector_f(SP_pix.nx);
  (*SP_lam).alocldo_flag=1;

  if(DEBUG) printf(" Aqui termine aloc\n");

  for(i=0;i<SP_pix.nx;i++) {
    if(DEBUG) printf(" Voy al %d\n",i);
    ldo1=(*SP_lam).ldomin+(i-0.5)*(*SP_lam).deltaldo;
    ldo2=(*SP_lam).ldomin+(i+0.5)*(*SP_lam).deltaldo;
    if(DEBUG) printf(" Despues ldo\n");
    pixreal1=pr_ldo2pix(ldo1,DP)-1;
    pixreal2=pr_ldo2pix(ldo2,DP)-1;
    if(DP.B>0) SWAP(pixreal1,pixreal2);
    if(DEBUG) printf(" antes ldo[]\n");
    (*SP_lam).ldo[i]=(ldo1+ldo2)/2.;
    if(DEBUG) printf(" Antes vista\n");
    if(DEBUG) printf(" pixreal1 %f real2 %f\n",pixreal1,pixreal2);
    (*SP_lam).spec[i]=sumspecpix_vista(SP_pix.spec,SP_pix.nx,pixreal1,pixreal2) * pr_dpdl((ldo1+ldo2)/2.,DP);
    if(DP.B>0) (*SP_lam).spec[i]=-(*SP_lam).spec[i];
    if(DEBUG) printf(" Factoe %f %f \n",sumspecpix_vista(SP_pix.spec,SP_pix.nx,pixreal1,pixreal2), pr_dpdl((ldo1+ldo2)/2.,DP));
    if(DEBUG) printf(" Final ldo %f  spec_ldo %f interp %f \n", (*SP_lam).ldo[i],(*SP_lam).spec[i],(pixreal2-pixreal1)*intspecpix_lin(SP_pix.spec,SP_pix.nx,(pixreal1+pixreal2)/2.,0));
  }
  if(DEBUG) printf(" Sali de SP_Pix\n");
}

void SP_ldo2pix(struct spectrum *SP_pix, struct spectrum SP_lam, struct disper_prism DP) {
  
  /* Paso de un espectro por lambda a un espectro por pixel */
  /* De paso alocatea el espectro de lambda */

  int i;
  float lambdareal1,lambdareal2;
  float pix1,pix2;
  float pixreal1,pixreal2;
  float pixmin,pixmax;

  pixmin=pr_ldo2pix(SP_lam.ldomin,DP);
  pixmax=pr_ldo2pix(SP_lam.ldomin+SP_lam.nx*SP_lam.deltaldo,DP);
  (*SP_pix).nx=(int)(pixmax-pixmin);
  (*SP_pix).npixels=(int)(pixmax-pixmin);
  if((*SP_pix).aloc_flag==1) free((*SP_pix).spec);
  (*SP_pix).spec=vector_f((*SP_pix).nx);
  (*SP_pix).aloc_flag=1;
  (*SP_pix).ldomin=1;
  (*SP_pix).deltaldo=1;
  if((*SP_pix).alocldo_flag==1) free((*SP_pix).ldo);
  (*SP_pix).alocldo_flag=1;


  for(i=0;i<(*SP_pix).nx;i++) {
    pix1=(float)(pixmin+i-0.5);
    pix2=(float)(pixmin+i+0.5);
    lambdareal1=pr_pix2ldo(pix1,DP);
    lambdareal2=pr_pix2ldo(pix2,DP);
    pixreal1=(lambdareal1-SP_lam.ldomin)/SP_lam.deltaldo;
    pixreal2=(lambdareal2-SP_lam.ldomin)/SP_lam.deltaldo;
    (*SP_pix).ldo[i]=i+1;
    (*SP_pix).spec[i]=sumspecpix_vista(SP_lam.spec,SP_lam.nx,pixreal1,pixreal2) * pr_dldp((pix1+pix2)/2.,DP);
/*     printf(" i %d lamdna %f ldo %f   PI %f LA %f\n",i,lambdareal,ldos[i],(*SP_pix).spec[i],SP_lam.spec[i]); */
  }
}
