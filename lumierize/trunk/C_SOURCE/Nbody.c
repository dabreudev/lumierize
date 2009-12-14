#include "modulos.h"

#define N 1000

struct xyzm {
  float x;
  float y;
  float z;
  float m;
};
struct vel {
  float vx;
  float vy;
  float vz;
};
struct  body {
  struct xyzm pos;
  struct vel velo;
  float pgas;
  float sfr;
  int nmerge;
};



void Force( struct xyzm xyzm[N], int k, float force[3], int nc);

void Init(struct body *galaxy,struct body *galaxy2);
int main()
{
  
  struct body galaxy[N],galaxy2[N];
  struct xyzm k1[N],k2[N],k3[N],k4[N],dum[N];
  struct vel v1[N],v2[N],v3[N],v4[N];
  float t,h,tmax;
  float force[3];
  int nstep,i,nc;
  printf(" Maximum time: ");
  scanf("%f",&tmax);

  printf("%f\n",tmax);
  printf(" Number of steps: ");
  scanf("%d",&nstep);
 
  h=tmax/nstep;
  cpgopen("?");
  Init(galaxy,galaxy2);
//  printf("Ha salifo bien\n");
  for(t=0;t<tmax;t+=h){
    printf("Tiempo %f\n",t);

    for(i=0;i<N;i++)  {
//      printf(" Galaxy     %d >> %f %f %f\n",i,galaxy[i].pos.x,galaxy[i].pos.y,galaxy[i].pos.m);
//      printf(" Galaxy vel %d >>  %f %f\n",i,galaxy[i].velo.vx,galaxy[i].velo.vy);
    cpgpt1(galaxy[i].pos.x,galaxy[i].pos.y,4);
  
    }
//    printf("Pause\n");
    
//    scanf("%f",&nul);
      
   nc=0;
    for(i=0;i<N;i++){
//      printf("POR AQUI PASA\n");
      if(galaxy[i].nmerge==i){
	dum[i].x=galaxy[i].pos.x;
	dum[i].y=galaxy[i].pos.y;
	dum[i].z=galaxy[i].pos.z;
	dum[i].m=galaxy[i].pos.m;
//	printf("dum %f %f %f",dum[i].x,dum[i].y,dum[i].z);
	nc++;
      }
    }
    for(i=0;i<N;i++){
      Force(dum,i,force,nc);
//      printf("Fuerza : %f %f %f\n",force[0],force[1],force[2]);
      k1[i].x=h*galaxy[i].velo.vx;
      k1[i].y=h*galaxy[i].velo.vy;
      k1[i].z=h*galaxy[i].velo.vz;
      v1[i].vx=h*force[0]/galaxy[i].pos.m;
      v1[i].vy=h*force[1]/galaxy[i].pos.m;
      v1[i].vz=h*force[2]/galaxy[i].pos.m;
    }
    nc=0;
    for(i=0;i<N;i++){
      if(galaxy[i].nmerge==i){
	dum[i].x=galaxy[i].pos.x+k1[i].x/2;
	dum[i].y=galaxy[i].pos.y+k1[i].y/2;
	dum[i].z=galaxy[i].pos.z+k1[i].z/2;
	dum[i].m=galaxy[i].pos.m;       
	nc++;
      }
    }

    for(i=0;i<N;i++){
      Force(dum,i,force,nc);
//      printf("Fuerza : %f %f %f\n",force[0],force[1],force[2]);
      k2[i].x=h*(galaxy[i].velo.vx+v1[i].vx/2);
      k2[i].y=h*(galaxy[i].velo.vy+v1[i].vy/2);
      k2[i].z=h*(galaxy[i].velo.vz+v1[i].vz/2);
      v2[i].vx=h*force[0]/galaxy[i].pos.m;
      v2[i].vy=h*force[1]/galaxy[i].pos.m;
      v2[i].vz=h*force[2]/galaxy[i].pos.m;
    }
    nc=0;
    for(i=0;i<N;i++){
      if(galaxy[i].nmerge==i){
	dum[i].x=galaxy[i].pos.x+k2[i].x/2;
	dum[i].y=galaxy[i].pos.y+k2[i].y/2;
	dum[i].z=galaxy[i].pos.z+k2[i].z/2;
	dum[i].m=galaxy[i].pos.m;
	nc++;
      }
    }
    for(i=0;i<N;i++){
      Force(dum,i,force,nc);
//      printf("Fuerza : %f %f %f\n",force[0],force[1],force[2]);
      k3[i].x=h*(galaxy[i].velo.vx+v2[i].vx/2);
      k3[i].y=h*(galaxy[i].velo.vy+v2[i].vy/2);
      k3[i].z=h*(galaxy[i].velo.vz+v2[i].vz/2);
      v3[i].vx=h*force[0]/galaxy[i].pos.m;
      v3[i].vy=h*force[1]/galaxy[i].pos.m;
      v3[i].vz=h*force[2]/galaxy[i].pos.m;
    }

    nc=0;
    for(i=0;i<N;i++){
      if(galaxy[i].nmerge==i){
	dum[i].x=galaxy[i].pos.x+k3[i].x/2;
	dum[i].y=galaxy[i].pos.y+k3[i].y/2;
	dum[i].z=galaxy[i].pos.z+k3[i].z/2;
	dum[i].m=galaxy[i].pos.m;
	nc++;
      }
    } 
    for(i=0;i<N;i++){
      Force(dum,i,force,nc);
//      printf("Fuerza : %f %f %f\n",force[0],force[1],force[2]);
      k4[i].x=h*(galaxy[i].velo.vx+v3[i].vx);
      k4[i].y=h*(galaxy[i].velo.vy+v3[i].vy);
      k4[i].z=h*(galaxy[i].velo.vz+v3[i].vz);
      v4[i].vx=h*force[0]/galaxy[i].pos.m;
      v4[i].vy=h*force[1]/galaxy[i].pos.m;
      v4[i].vz=h*force[2]/galaxy[i].pos.m;
    }
    for(i=0;i<N;i++){
      galaxy2[i].velo.vx=galaxy[i].velo.vx+v1[i].vx/6+v2[i].vx/3+v3[i].vx/3+v4[i].vx/6;
      galaxy2[i].velo.vy=galaxy[i].velo.vy+v1[i].vy/6+v2[i].vy/3+v3[i].vy/3+v4[i].vy/6;
      galaxy2[i].velo.vz=galaxy[i].velo.vz+v1[i].vz/6+v2[i].vz/3+v3[i].vz/3+v4[i].vz/6;
      galaxy2[i].pos.x=galaxy[i].pos.x+k1[i].x/6+k2[i].x/3+k3[i].x/3+k4[i].x/6;
      galaxy2[i].pos.y=galaxy[i].pos.y+k1[i].y/6+k2[i].y/3+k3[i].y/3+k4[i].y/6;
      galaxy2[i].pos.z=galaxy[i].pos.z+k1[i].z/6+k2[i].z/3+k3[i].z/3+k4[i].z/6;
//      printf(" 2322 Galaxy %d\n  %f %f %f",i,galaxy2[i].pos.x,galaxy2[i].pos.y,galaxy2[i].pos.m);
//      printf(" Galaxy vel%d  %f %f\n",i,galaxy2[i].velo.vx,galaxy2[i].velo.vy);
    }
    cpgsci(0);
    cpgpt1(galaxy[i].pos.x,galaxy[i].pos.y,4);
    cpgsci(1);


//    printf("Ha acabdodo \n");    

    memcpy(galaxy,galaxy2,N*sizeof(struct body));
  }


  return 0;
}


void Force( struct xyzm xyzm[N], int k, float force[3], int nc) {
  float G=39.478418;     // Para el sistema solar. Unidades: UA, año, Msolar.
  float xx,yy,zz,r,factor;
  int i;
  force[0]=0;force[1]=0;force[2]=0;
//      printf("YES %d\n",nc);

  for(i=0;i<nc;i++) {
//    printf("YES\n");
    if(i!=k){
      xx=xyzm[i].x-xyzm[k].x;
      yy=xyzm[i].y-xyzm[k].y;
      zz=xyzm[i].z-xyzm[k].z;
//      printf("xx %f yy %f zz %f\n",xx,yy,zz);
      r=sqrt(xx*xx+yy*yy+zz*zz);
      factor=xyzm[i].m*xyzm[k].m/(r*r*r);
      force[0]=factor*xx+force[0];  
      force[1]=factor*yy+force[1];   
      force[2]=factor*zz+force[2];
//      printf("%f <<<<<<<<\n",force[0]);
//      printf("%f <<<<<<<<\n",force[1]);
//      printf("%f <<<<<<<<\n",force[2]);
   }
  }
  force[0]=G*force[0]; 
  force[1]=G*force[1];   
  force[2]=G*force[2];
}



void Init(struct body *galaxy,struct body *galaxy2) 
{
  int i;
  int op;
  float m,dm,vm,dvm,l;
  printf(" 1 Ramdomly distribution\n 2 Specify distribution\n");
  printf("\n Input option? ");
  scanf("%d",&op); 
  srandom((unsigned int)time(NULL)/2); 

  if(op==1) {
    printf(" Maximum lenght: ");
    scanf("%f",&l);
    printf(" Mean mass: ");
    scanf("%f",&m);
    printf(" Stdv. mass: ");	
    scanf("%f",&dm);
    printf(" Mean vel; ");
    scanf("%f",&vm);
    printf(" Stdv. vel: ");	
    scanf("%f",&dvm);
    cpgsvp(0.,1.,0.,1.);
    cpgswin(-l*4,l*4,-l*4,l*4);
    cpgsch(1.);
//    printf("random: %f\n",random()/2147483647.);
    for(i=0;i<N;i++) 
      {
	galaxy[i].pos.x=random()/2147483647.*2*l-l;
	galaxy[i].pos.y=random()/2147483647.*2*l-l;
	galaxy[i].pos.z=0;
//	galaxy[i].pos.z=random()/2147483647.*2*l-l;
	galaxy[i].pos.m=random()/2147483647.*dm+m;
	galaxy[i].nmerge=i;
	galaxy2[i].pos.m=galaxy[i].pos.m;
	galaxy2[i].nmerge=galaxy[i].nmerge;
	galaxy[i].velo.vx=random()/2147483647.*vm-vm/2;
	galaxy[i].velo.vy=random()/2147483647.*vm-vm/2;
	galaxy[i].velo.vz=0;
//	galaxy[i].velo.vz=random()/2147483647.*vm-vm/2;
	
      }
  }
  if(op==2) {
    for(i=0;i<N;i++) {
      printf(" Particle %d\n",i);
      printf(" X: ");
      scanf("%f",&(galaxy[i].pos.x));
      printf(" Y: ");
      scanf("%f",&(galaxy[i].pos.y));
      printf(" Z: ");
      scanf("%f",&(galaxy[i].pos.z));
      printf(" Vel X: ");
      scanf("%f",&(galaxy[i].velo.vx));
      printf(" Vel Y: ");
      scanf("%f",&(galaxy[i].velo.vy));
      printf(" Vel Z: ");
      scanf("%f",&(galaxy[i].velo.vz));
      printf(" Mass: ");
      scanf("%f",&(galaxy[i].pos.m));
      galaxy[i].nmerge=i;
      galaxy2[i].pos.m=galaxy[i].pos.m;
      galaxy2[i].nmerge=galaxy[i].nmerge;
      
    }
    printf(" Maximum lenght: ");
    scanf("%f",&l);
    cpgsvp(0.,1.,0.,1.);
    cpgswin(-l*2,l*2,-l*2,l*2);
    cpgsch(1.);
    
    
  }
  
}
