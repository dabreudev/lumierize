 
#include "modulos.h" 

void Gausfilter2D(float *ima_input, float *ima_output,int nx, int ny, float sigma) 
{ 
  /* Aplica un filtro gaussiano a una imagen bidimensional. 
     La imagen de entrada es ima_input tiene nx*ny pixels, y la anchura de la gaussiana es sigma pixels 
     La imagen de salida es ima_output de nx*ny pixels y 
     DEBE ESTAR DIMENSIONADA COMO nx*ny 
     
     */ 
  
  int i,j; 
  int k,l; 
  float expo;
  float pi=4.*atan(1.);
/*   float tr[6]; */
/*   float min,max; */
  if(sigma==0) {
/*     printf("Ha pasado por aui %f\n",pi); */
    for(i=0;i<nx*ny;i++) { 
      ima_output[i]=ima_input[i];
    }
    return;
  }
  /*
  min=max=0;
  cpgpanl(2,2);
  cpgswin(0,nx,0,ny);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  printf("aqui no\n");
  cpggray(ima_input,nx,ny,1,nx-1,1,ny-1,1e-21,min,tr);
  */
  for(i=0;i<nx;i++) { 
    for(j=0;j<ny;j++) { 
      *(ima_output+i+j*nx)=0;
      for(k=0;k<nx;k++) { 
	for(l=0;l<ny;l++) { 
/* 	  printf("por qaqui va bien %d %d %d %d\n",i,j,k,l); */
	  expo=-((i-k)*(i-k)+(j-l)*(j-l))/2./sigma/sigma; 
/* 	  printf("expo %f  %f %f %f %f\n",expo,*(ima_input+k+l*nx),exp(expo),*(ima_input+k+l*nx)*exp(expo),*(ima_output+i+j*nx)); */
	  *(ima_output+i+j*nx)+=*(ima_input+k+l*nx)*exp(expo); 
/* 	  printf("por qaqui va bien %f %f\n",exp(expo),expo); */
	} 
      } 
/*       printf("El valor de salida es %f \n",*(ima_output+i*nx+j)); */
      (*(ima_output+i+j*nx))*=1./2./pi/sigma/sigma; 
/*       if(*(ima_output+i*nx+j)< min) min=*(ima_output+i*nx+j); */
/*       if(*(ima_output+i*nx+j)> max) max=*(ima_output+i*nx+j); */
/*       printf("El valor de salida es %f \n",*(ima_output+i*nx+j)); */
   } 
  } 
   /*cpgpanl(2,1);
  cpgswin(0,nx,0,ny);
  cpgbox("BCTNS",0,0,"BCTNS",0,0);
  tr[0]=0;tr[1]=1;tr[2]=0;tr[3]=0;tr[4]=0;tr[5]=1;
  printf("min  %f max %f\n",min,max);
  cpggray(ima_output,nx,ny,1,nx-1,1,ny-1,max,min,tr);
  */
} 

