#include "modulos.h"



int main(int argc, char **argv) {


  char filename[100];
  char c1,c2,c3;
  int flag=0;
  int nchar,i;
  int i1,i2;
  char comando[200];
  
  FILE *f;
  FILE *fo;
  if(argc < 2) {
    printf(" Use: converpar parfile\n");
    exit(1);
  }

  strcpy(filename,argv[1]);

  f=fopen(filename,"r");
  fo=fopen("fichero_temporal_para_convertpar","w");
 
  nchar=0;

  do {
    c1=fgetc(f);
    if(c1=='E') {
      c2=getc(f);
      if(c2=='N') {
	c3=getc(f);
	if(c3=='D') {
	  flag=1;
	}
	else {
	  fputc(c1,fo);
	  fputc(c2,fo);
	  fputc(c3,fo);
	  nchar++;
	  nchar++;
	  nchar++;
	}
	
      }
      else {
	fputc(c1,fo);
 	fputc(c2,fo);
	nchar++;
	nchar++;
      }
    }
    else {
      fputc(c1,fo);
      nchar++;
    }
    

 }while(!flag);

  
  printf("There are %d 2880 blocks + trash",((int)(nchar/2880)));

  printf(" Filling from %d to %d with blanks\n",nchar,((int)(nchar/80)+1)*80);

  i1=nchar;i2=((int)(nchar/80)+1)*80-1;
  
  for(i=i1; i< i2;i++) {
    fputc(' ',fo);
    nchar++;
  }
  fputc('\n',fo);nchar++;

  printf(" Filling from line %d to line %d with COMMENT lines\n",(int)(nchar/80),((int)(nchar/2880)+1)*36);

  i1=(int)(nchar/80);i2=((int)(nchar/2880)+1)*36;

  for(i=i1; i< i2-1;i++) {
    fputc('C',fo);fputc('O',fo);fputc('M',fo);fputc('M',fo);fputc('E',fo);fputc('N',fo);fputc('T',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc('\n',fo);
    nchar+=80;
  }

  
  
  
  fputc('E',fo);fputc('N',fo);fputc('D',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc(' ',fo);fputc('\n',fo);
  
  
  fclose(f);
  fclose(fo);

  sprintf(comando,"mv %s %s.old\n",filename,filename);
  system(comando);
  sprintf(comando,"mv fichero_temporal_para_convertpar %s\n",filename);
  system(comando);
  sprintf(comando,"rm -f  fichero_temporal_para_convertpar\n");
  system(comando);

  return 0;
}
