#include "modulos.h"

float area_in(float x0,float y0,float r,float xc,float yc,float a);
float int_circ(float r,float x1,float x2);


float circ_aper(float *ima, int nx, int ny,float x,float y,
		float r,int *npix){
  /*
    x:     Centro de la apertura. CUIDADO. Este centro esta dado en  
    y:     unidades de la matriz. Normalmente, para una imagen, no corres
           ponde con los pixeles. Si quiero hallar la fotometria de apertura
           en el pixel xp, yp, tendre que hacer una llamada a la funcion como
           circ_aper(           xp-1,yp-1,..      )
  */
  register int i,j;
  int i1,i2,j1,j2;
  float xm,ym;
  float s=0;

  if(0.5>r){
    i=(int)floor(x);
    j=(int)floor(y);
    printf("Atencion, la apertura es menor que el pixel\n");
    return M_PI*r*r*ima[j*nx+i];
  }

  i1=(int)(x-r-2);
  i2=(int)(x+r+2);
  j1=(int)(y-r-2);
  j2=(int)(y+r+2);
  if(i1<0) i1=0;
  if(j1<0) j1=0;
  if(i2>nx-1) i2=nx-1;
  if(j2>ny-1) j2=ny-1;
  for(j=j1;j<=j2;j++) {
    for(i=i1;i<=i2;i++) {
      xm=i+0.5;
      ym=j+0.5;
      s+=ima[j*nx+i]*area_in(x,y,r,xm,ym,1);
/*       printf("aper %d %d : ima %f fac %f s %f \n",i,j,ima[j*nx+i],area_in(x,y,r,xm,ym,1),s); */
/*      area+=area_in(x,y,r,xm,ym,1); */
    }
  }
  *npix=(int)(M_PI*r*r);
  return s;
}


float  circ_aper_err(float *ima, float *errima, int nx, int ny, 
		     float x,float y,float r,int *npix, float *errflux){
  
  /*
    x:     Centro de la apertura. CUIDADO. Este centro esta dado en  
    y:     unidades de la matriz. Normalmente, para una imagen, no corres
    ponde con los pixeles. Si quiero hallar la fotometria de apertura
    en el pixel xp, yp, tendre que hacer una llamada a la funcion como
    circ_aper(           xp-1,yp-1,..      )
    esta subrutina calcula ademas los errores en la fotometria dada 
    una imagen de errores
  */

  register int i,j;
  int i1,i2,j1,j2;
  float xm,ym;
  float s=0;
  *errflux=0;

  if(0.5>r){
    i=(int)floor(x);
    j=(int)floor(y);
    printf("Atencion, la apertura es menor que el pixel\n");
    return M_PI*r*r*ima[j*nx+i];
  }
  i1=(int)(x-r-2);
  i2=(int)(x+r+2);
  j1=(int)(y-r-2);
  j2=(int)(y+r+2);
  if(i1<0) i1=0;
  if(j1<0) j1=0;
  if(i2>nx-1) i2=nx-1;
  if(j2>ny-1) j2=ny-1;
  for(j=j1;j<=j2;j++) {
    for(i=i1;i<=i2;i++) {
      float factor;
      xm=i+0.5;
      ym=j+0.5;
      factor=area_in(x,y,r,xm,ym,1);
      s+=ima[j*nx+i]*factor;
/*       printf("aper %d %d x %f y %f : ima %f fac %f s %f \n",i,j,ima[j*nx+i],factor,s); */
      (*errflux)+=errima[j*nx+i]*errima[j*nx+i]*factor*factor;
/*      area+=area_in(x,y,r,xm,ym,1);   */
    }
  }
  (*errflux)=sqrt(*errflux);
  *npix=(int)(M_PI*r*r);
/*  *npix=(int)area;  */
  return s;
}


float area_in(float x0,float y0,float r,float xc,float yc,float a){
  /* (r,x0,y0) radi y coor del centro del circulo
     (a,xc,yc) lado y coor del centro del cuadrado
     devuelve el area dentro del cuadrado si r>a/2
  */
  float lad=a/2;
  if(lad>r){
    printf("Atenciorr, la apertura es menor que el pixel\n");
    return -1;
  }
  else{
    float axc=x0+fabs(x0-xc);
    float ayc=y0+fabs(y0-yc);
    float y1=ayc-lad;
    float y2=ayc+lad;
    float x1=axc-lad;
    float x2=axc+lad;
    float a=(x2-x1)*(x2-x1);
    float b=2*(x1-x0)*(x2-x1);
    float c1=(x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)-r*r;
    float c2=(x1-x0)*(x1-x0)+(y2-y0)*(y2-y0)-r*r;
    float delta1=b*b-4*a*c1;
    float delta2=b*b-4*a*c2;
    float t1,t2,u1,u2;
    if(delta1<0){ /* esta por fuera y por arriba */
      return 0;
    }
    else{
      t1=(-b+sqrt(delta1))/(2*a);
      t2=(-b-sqrt(delta1))/(2*a);
      if(t1>=1 && delta2>=0){ /* dentro o borde */
	u1=(-b+sqrt(delta2))/(2*a);
	u2=(-b-sqrt(delta2))/(2*a);
	if(u1>=1){ /* dentro */
	  return a*a;
	}
	else{
	  if(u1>0){
	    if(u2>0){ /* borde dos cortes arriba */
	      float cut1=x1+u2*(x2-x1);
	      float cut2=x1+u1*(x2-x1);
 	      return (a*a*(int_circ(r,x1-x0,cut1-x0)+(cut2-cut1)*(y2-y0)+int_circ(r,cut2-x0,x2-x0)-(y1-y0)*(x2-x1))); 
	    }
	    else{ /* borde un corte arriba */
	      float cut1=x1+u1*(x2-x1);
	      return (a*a*((cut1-x1)*(y2-y1)+int_circ(r,cut1-x0,x2-x0)-(y1-y0)*(x2-cut1)));
	    }
	  }
	  else{ /* Extremos borde la caja */
	    return a*a*int_circ(r,x1-x0,x2-x0)-(y1-y0)*(x2-x1);
	  }
	}
      }
      if(t1>=1 && delta2<0) /* Extremos borde la caja */
	return a*a*(int_circ(r,x1-x0,x2-x0)-(y1-y0)*(x2-x1));
      if(t1<=0)  /* fuera */
	return 0;
      if(t1>=0 && t1<=1 && delta2<0){
	if(t2>=0 && t2<=1){  /* borde. Dos cortes abajo Ninguno arriba*/
	  float cut1=x1+t2*(x2-x1);
	  float cut2=x1+t1*(x2-x1);
	  return a*a*(int_circ(r,cut1-x0,cut2-x0)-(y1-y0)*(cut2-cut1));
	}
	else{ /* borde. Un corte abajo. Ninguno arriba */
	  float cut1=x1+t1*(x2-x1);
	  return a*a*(int_circ(r,x1-x0,cut1-x0)-(y1-y0)*(cut1-x1));
	}
      }
      if(t1>=0 && t1<=1 && delta2>0){
	u1=(-b+sqrt(delta2))/(2*a);
	u2=(-b-sqrt(delta2))/(2*a);
	if(t2>=0 && t2<=1){
	  if(u1>=0 && u1<=1){
	    if(u2>=0 && u2<=1){ /*Borde. Dos cortes abajo y dos arriba */
	      float cut1=x1+t2*(x2-x1);
	      float cut2=x1+t1*(x2-x1);
	      float up_cut1=x1+u2*(x2-x1);
	      float up_cut2=x1+u1*(x2-x1);
	      return a*a*(int_circ(r,cut1-x0,up_cut1-x0)
			  +(y1-y0)*(up_cut2-up_cut1)+
			  int_circ(r,up_cut2-x0,cut1-x0)-(y1-y0)*(cut2-cut1));
	    }
	    else{ /*Borde. Dos cortes abajo y uno arriba */
	      printf("Never reaches this point 1\n");
	      return 0;
	    }
	  }
	  else{ /* borde. Dos cortes abajo */
	    float cut1=x1+u2*(x2-x1);
	    float cut2=x1+u1*(x2-x1);
	    return a*a*(int_circ(r,cut1-x0,cut2-x0)-(y1-y0)*(x2-x1));
	  }
	}
	else{
	  if(u1>=0 && u1<=1){ /* Borde. Un corte abajo y uno arriba */
	    return a*a*(int_circ(r,y1-y0,y2-y0)-(x1-x0)*(y2-y1));
	  }
	  else{ /* Borde. Un corte abajo */
	    float cut1=x1+t1*(x2-x1);
	    return a*a*(int_circ(r,x1-x0,cut1-x0)-(y1-y0)*(cut1-x1));
	  }
	}
      }
    }
    printf("Never reaches this point 2\n");
    return 0;
  }
}

float int_circ(float r,float x1,float x2){
  float a1=x1/r;
  float a2=x2/r;
  return r*r*(asin(a2)+a2*sqrt(1-a2*a2)-asin(a1)-a1*sqrt(1-a1*a1))/2;
}


float  elip_aper(float *ima, int nx, int ny, float x,float y,
		 float a, float b, float pa,int *npix){
  int i,j;
  int i1,i2,j1,j2;
  float flux=0;
  (*npix)=0;

/*   //El angulo de posicion tiene que venir dado en radianes!! */
  
  i1=(int)(x-b);
  i2=(int)(x+b);
  j1=(int)(y-b);
  j2=(int)(y+b);
  if(i1<0) i1=0;
  if(j1<0) j1=0;
  if(i2>nx-1) i2=nx-1;
  if(j2>ny-1) j2=ny-1;
  for(i=i1;i<=i2;i++)  {
    for(j=j1;j<=j2;j++) {
      if( (((i-x)*cos(pa)+(j-y)*sin(pa))*((i-x)*cos(pa)+(j-y)*sin(pa))/a/a+(-(i-x)*sin(pa)+(j-y)*cos(pa))*(-(i-x)*sin(pa)+(j-y)*cos(pa))/b/b)<=1 ) {
	(*npix)++;
	flux+=ima[i+nx*j];
      }
    }
  }
  
  return flux;
}



/* Rutinas antiguas */


float  circ_aper_bak(float *ima, int nx, int ny, 
		     float x,float y,float r,int *npix)
{
  
  /*
    x:     Centro de la apertura. CUIDADO. Este centro esta dado en  
    y:     unidades de la matriz. Normalmente, para una imagen, no corres
           ponde con los pixeles. Si quiero hallar la fotometria de apertura
           en el pixel xp, yp, tendre que hacer una llamada a la funcion como
           circ_aper(           xp-1,yp-1,..      )
  */

  /* float xcen,ycen; */

  int i,j;
  int i1,i2,j1,j2;
  float flux=0;
  /* printf(" Entra en suvr\n"); */
  (*npix)=0;
  /* printf(" ASIGNA  NX %d NY %d\n",nx,ny); */
  
  i1=(int)(x-r-2);
  i2=(int)(x+r+2);
  j1=(int)(y-r-2);
  j2=(int)(y+r+2);
  if(i1<0) i1=0;
  if(j1<0) j1=0;
  if(i2>nx-1) i2=nx-1;
  if(j2>ny-1) j2=ny-1;
  /*  printf(" El radio es %f y el centro %f,%f\n",r,x,y); */
  /*  printf(" Los limites son : i: %d %d   j:%d %d \n",i1,i2,j1,j2); */
  /*  xcen=0;ycen=0; */
  for(i=i1;i<=i2;i++) {
    for(j=j1;j<=j2;j++) {
      /* printf(" i %d j %d",i,j); */
      if (((float)i-x)*((float)i-x)+((float)j-y)*((float)j-y)<=r*r) {
	/* xcen+=i*ima[i+nx*j]; */
	/* ycen+=j*ima[i+nx*j]; */
	(*npix)++;
	flux+=ima[i+nx*j];
	/* printf(" npix %d flux %f  ima %f\n",*npix,flux,ima[i+nx*j]); */
      }
    }
  }
  /* xcen/=flux; */
  /* ycen/=flux; */
  /* printf("EN %f %f CAL %f %f  DIFF %f %f\n",x,y,xcen,ycen,x-xcen,y-ycen); */
  return(flux);
}



float  circ_aper_err_bak(float *ima, float *errima, int nx, int ny, float x,float y,float r,int *npix, float *errflux)
{
  
  /*
    x:     Centro de la apertura. CUIDADO. Este centro esta dado en  
    y:     unidades de la matriz. Normalmente, para una imagen, no corres
           ponde con los pixeles. Si quiero hallar la fotometria de apertura
           en el pixel xp, yp, tendre que hacer una llamada a la funcion como
           circ_aper(           xp-1,yp-1,..      )
    esta subrutina calcula ademas los errores en la fotometria dada una imagen de errores
  */

/*   //float xcen,ycen; */

  int i,j;
  int i1,i2,j1,j2;
  float flux=0;

  *errflux=0.;
  (*npix)=0;
  i1=(int)(x-r-2);
  i2=(int)(x+r+2);
  j1=(int)(y-r-2);
  j2=(int)(y+r+2);
  if(i1<0) i1=0;
  if(j1<0) j1=0;
  if(i2>nx-1) i2=nx-1;
  if(j2>ny-1) j2=ny-1;
  for(i=i1;i<=i2;i++) {
    for(j=j1;j<=j2;j++) {
      if (((float)i-x)*((float)i-x)+((float)j-y)*((float)j-y)<=r*r) {
	(*npix)++;
	flux+=ima[i+nx*j];
	(*errflux)+=errima[i+nx*j]*errima[i+nx*j];
      }
    }
  }
  (*errflux)=sqrt(*errflux);
  return(flux);
}
