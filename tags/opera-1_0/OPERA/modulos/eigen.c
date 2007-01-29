#include "modulos.h"

#define MAXITE 999

float modulo(int n,float *u);


void eigen(float *m,int l,float *eigenval,float *eigenvec)
{
  /*
    m:        square matrix of dimensions lxl
    l:        dimension of the matrix
    eigenval: output vector with eigenvalues. Dimension l
    eigenvec: output vector with eigenvectors. Each eigenvector is 
    a column.
  */

/*   int i,j,k,n,it,ita[l]; */
/*   float a[l][l],e,dif[l],per[l],dist; */
/*   float tr[MAXITE+1][l][l],pmo[l],dum,vat; */
  int i,j,k,it,*ita;
  float *a,e,*dif,*per,dist;
  float *tr,*pmo,dum,vat;
/*   int seed; */
  if((ita=malloc(l*sizeof(int)))==NULL)      printf("I cannot dimension ita of %d elements \n",l);
  if((a=malloc(l*l*sizeof(float)))==NULL)    printf("I cannot dimension a of %d elements \n",l*l);
  if((dif=malloc(l*sizeof(float)))==NULL)    printf("I cannot dimension dif of %d elements \n",l);
  if((per=malloc(l*sizeof(float)))==NULL)    printf("I cannot dimension per of %d elements \n",l);
  if((pmo=malloc(l*sizeof(float)))==NULL)    printf("I cannot dimension pmo of %d elements \n",l);
  if((tr=malloc((MAXITE+1)*l*l*sizeof(float)))==NULL) printf("I cannot dimension tr of %d elements \n",l*l*l);


/*   for(i=0;i<l;i++) { */
/*     printf(" %f %f %f %f %f\n",m[0+l*i],m[1+l*i],m[2+l*i],m[3+l*i],m[4+l*i]); */
/*     printf(" %f %f %f %f %f\n",y[0+5*i],y[1+5*i],y[2+5*i],y[3+5*i],y[4+5*i]); */
/*   } */
  srandom((unsigned int)time(NULL)/2);
/*   //printf(" l %d\n",l); */
  e=1.e-9;
/*   //seed=-2; */
/*   //printf(" Este es el principi\n"); */
/*   //exit(1); */
  for(j=0;j<l;j++) {
    for(k=0;k<l;k++) {
/*       //a[j][k]=m[j+l*k]; */
      a[j+l*k]=m[j+l*k];
    }
  }
  vat=0.;
  for(i=0;i<l;i++) {            /* //eigenvalues loop */
/*     //printf(" Eigenval %d\n",i); */
    eigenval[i]=0.;
/*     //printf(" Matriz \n"); */
    for(j=0;j<l;j++) {
      for(k=0;k<l;k++) {
/* 	//if(!(i==0)) a[j][k]=a[j][k]-eigenval[i-1]*eigenvec[j+l*(i-1)]*eigenvec[k+l*(i-1)]; */
	if(!(i==0)) a[j+l*k]=a[j+l*k]-eigenval[i-1]*eigenvec[j+l*(i-1)]*eigenvec[k+l*(i-1)];
/* 	//printf(" %g ",a[j+l*k]); */
      }
/*       //printf("\n"); */
    }
    
    for(j=0;j<l;j++) {
      for(k=0;k<MAXITE;k++) {
/* 	//tr[k][j][i]=1+(2*random()/2147483647.-1.)/10.; */
	tr[k+(MAXITE+1)*j+(MAXITE+1)*l*i]=1+(2*random()/2147483647.-1.)/10.;
      }
      eigenvec[j+l*i]=0.;
    }
/*     //printf(" Aqui no\n"); */
    for(it=0;it<MAXITE;it++) {
      dum=0.;
      for(j=0;j<l;j++) {
	dum=0.;
	for(k=0;k<l;k++) {
/* 	  //dum=dum+a[j][k]*tr[it][k][i]; */
	  dum=dum+a[j+l*k]*tr[it+(MAXITE+1)*k+(MAXITE+1)*l*i];
	}
       
	eigenvec[j+l*i]=dum;
      }
      for(j=0;j<l;j++) {
	pmo[j]=eigenvec[j+l*i];
      }
      for(j=0;j<l;j++) {
/* 	//tr[it+1][j][i]=eigenvec[j+l*i]/modulo(l,pmo); */
	tr[it+1+(MAXITE+1)*j+(MAXITE+1)*l*i]=eigenvec[j+l*i]/modulo(l,pmo);
      }
      if(!(it==1)) {
	dist=0.;
	for(j=0;j<l;j++) {
/* 	  //dist=dist+(tr[it+1][j][i]-tr[it][j][i])*(tr[it+1][j][i]-tr[it][j][i]); */
	  dist=dist+(tr[it+1+(MAXITE+1)*j+(MAXITE+1)*l*i]-tr[it+(MAXITE+1)*j+(MAXITE+1)*l*i])*(tr[it+1+(MAXITE+1)*j+(MAXITE+1)*l*i]-tr[it+(MAXITE+1)*j+(MAXITE+1)*l*i]);
	}
	dist=sqrt(dist);
	if(dist<e) break;
      }
      
    }
/*     //printf(" Aqui tambpoco\n"); */
    ita[i]=it;
/*     //printf(" EIGENVEC  "); */
    for(j=0;j<l;j++) {
/*       //printf("  %f ",eigenvec[j+l*i]); */
      eigenvec[j+l*i]=eigenvec[j+l*i]/modulo(l,pmo);
    }
/*     //printf(" AQIo si \n"); */
/*     //printf("\nMODULO %f\n",modulo(l,pmo)); */
    
    eigenval[i]=modulo(l,pmo);
    vat=vat+eigenval[i];
    printf(" Eigenval %d: %f\n",i,eigenval[i]);
  }
/*   for(i=0;i<l;i++) { */
/*     per[i]=eigenval[i]*100/vat; */
/*     printf(" autoval %d en iter %d: %g ( %f)\n",i,ita[i],eigenval[i],per[i]); */
/*     printf(" autovec %d: ",i); */
/*     for(j=0;j<l;j++) { */
/*       printf(" %f ",eigenvec[j+l*i]); */
/*     } */
/*     printf("\n"); */
/*   } */
}


float modulo(int n,float *u) {
  int i;
  float dum;
  dum=0.;
  for(i=0;i<n;i++) {
    dum+=u[i]*u[i];
  }
  dum=sqrt(dum);
  return(dum);
}




/*        SUBROUTINE EIGEN(M,L,VA,VE) */
/* C A  Matriz cuadrada de datos */
/* C L  Dimension de la matriz */
/* C VA Vector de autovalores (VA(1),1er autovalor, VA(2), 2o autovalor, etc) */
/* C VE Matriz con los autovectores por columnas (VE(*,1), 1er autovector, etc) */
/*        INTEGER L,I,J,K,N,ITE,IT,ITA(L) */
/*        PARAMETER(ITE=199) !NUMERO MAXIMO DE ITERACIONES */
/*        REAL M(L,L),A(L,L),VA(L),VE(L,L),E,DIF(L),PER(L) */
/*        REAL MODULO,TR(ITE+1,L,L),PMO(L),DUM,VAT,RANRED */
/*        INTEGER SEED */
/*        EXTERNAL MODULO,RANRED */
/*        E=1.E-6 !LIMITE DE LA ITERACION */
/*        SEED=-2 */

/*        DO J=1,L */
/*        DO K=1,L */
/*        A(J,K)=M(J,K) */
/*           ENDDO */
/*        ENDDO */

/*        VAT=0.0 */
/*        DO I=1,L !BUCLE DE NUMERO DE AUTOVALORES */
/*        VA(I)=0.0 !INICIALIZACION DE LOS AUTOVALORES */

/*        DO J=1,L */
/*        DO K=1,L */
/*        IF (I.NE.1) A(J,K)=A(J,K)-VA(I-1)*VE(J,I-1)*VE(K,I-1) !REDEFIN. DE A  */
/*        ENDDO */
/*        ENDDO */

/*        DO J=1,L */
/*        DO K=1,ITE */
/*        TR(K,J,I)=1.+(2.*RANRED(SEED)-1.)/10. */
/*        ENDDO */
/*        VE(J,I)=0.        */
/*        ENDDO */

/*        DO 5 IT=1,ITE */
/*        DUM=0. */

/*        DO J=1,L */
/*        DUM=0. */
/*        DO K=1,L */
/*        DUM=DUM+A(J,K)*TR(IT,K,I) */
/*        ENDDO */
/*        VE(J,I)=DUM */
/*        ENDDO */

/*        DO J=1,L */
/*        PMO(J)=VE(J,I) */
/*        ENDDO */

/*        DO J=1,L */
/*        TR(IT+1,J,I)=VE(J,I)/MODULO(L,PMO) */
/*        ENDDO */

/*        IF (IT.NE.1) THEN */
/*        DIST=0. */
/*        DO J=1,L */
/*        DIST=DIST+(TR(IT+1,J,I)-TR(IT,J,I))**2. */
/*        ENDDO */
/*        DIST=SQRT(DIST) */
/*        IF (DIST.LE.E) GOTO 10 */
/*        ENDIF       */


/*  5     CONTINUE */

/*  10    CONTINUE */

/*  100   FORMAT(A10,1x,I1,A9,1x,I3,A6,F6.3,A2,F6.2,A2) */
/*  110   FORMAT(A9,3(A1,F6.3),A1) */
/*        ITA(I)=IT */
/*        DO J=1,L */
/*        VE(J,I)=VE(J,I)/MODULO(L,PMO) */
/*        ENDDO */

/*        VA(I)=MODULO(L,PMO) */
/*        VAT=VAT+VA(I) */
/*        ENDDO !FIN DE BUCLE EN AUTOVALORES */

       
       
/*        DO I=1,L */
/*        PER(I)=VA(I)*100./VAT */
/*        WRITE(*,100)  */
/*      > 'AUTOVALOR ',I,' EN ITER.',ITA(I),' VALE ',VA(I),' (',PER(I),'%)' */
/*        WRITE(*,110) 'AUTOVECT=','(',VE(1,I),',',VE(2,I),',',VE(3,I),')' */
/*        ENDDO */


/*        RETURN */
/*        END */

/*        REAL FUNCTION MODULO(N,U) */
/*        INTEGER N,I */
/*        REAL U(N),DUM */
/*        MODULO=0. */
/*        DO I=1,N */
/*        DUM=MODULO+U(I)*U(I) */
/*        MODULO=MODULO+U(I)*U(I) */
/*        ENDDO */
/*        MODULO=SQRT(MODULO) */
/*        RETURN */
/*        END */

