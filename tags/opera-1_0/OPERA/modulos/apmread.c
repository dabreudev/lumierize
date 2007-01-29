/*
 * Get data from the APM catalogues server. 
 *
 * This is a standalone program which makes a request of the
 * APM catalogues system based on command line arguments and returns
 * either a list or a postscript finding chart on standard output.
 *
 * VMS Notes:
 *    This program assumes Multinet TCP/IP support.
 *    This program expects to run as a foreign command.
 *    If standard output (SYS$OUTPUT) is redirected to a file (using
 *      ASSIGN or DEFINE) the resulting file will have control
 *      characters before each record.  
 *   
 * This version should compile under Ultrix and VMS/Multinet.
 *
 * Created: 1-February-1995 by T. McGlynn
 *          Goddard Space Flight Center
 *          Universities Space Research Association
 *          Code 668.1
 * Modified by Geraint Lewis and Mike Irwin to support APM online catalogues
 *          1-April-1996
 *
 */
/*#define APMCAT "131.111.68.247" - alternative form if name resolver is crap*/
/*#define APMCAT "www.aao.gov.au"*/
#define APMCAT "www.ast.cam.ac.uk"
#define PORT	80
#define ERROR  -1

#include "modulos.h"


/* Define return types of functions */
char	*mo_code(char *string);
char	*form_request(char *coords, char *survey, char *equinox, char *list, 
         char *box, char *numbers,char *img, char *email);
/* //char	*strchr(); */
int      getbuf(int s, char *rbuf, int qlen);



/* Request buffer */
char	req_str[2048];

char	*out_file = 0;

/* Is this an informational request */
int	info_req = 0;

int apmread(cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,classd,nstarmax,
	    gnum,gra,gdec,gmag,gtype,nlog)
     
     double  cra;            /* Search center J2000 right ascension in degrees */
     double  cdec;           /* Search center J2000 declination in degrees */
     double  dra;            /* Search half width in right ascension in degrees */
     double  ddec;           /* Search half-width in declination in degrees */
     double  drad;           /* Limiting separation in degrees (ignore if 0) */
     int     sysout;         /* Search coordinate system */
     double  eqout;          /* Search coordinate equinox */
     double  epout;          /* Proper motion epoch (0.0 for no proper motion) */
     double  mag1,mag2;      /* Limiting magnitudes (none if equal) */
     int     classd;         /* Desired object class (-1=all, 0=stars, 3=nonstars) */
     int     nstarmax;       /* Maximum number of stars to be returned */
     double  *gnum;          /* Array of APM numbers (returned) */
     double  *gra;           /* Array of right ascensions (returned) */
     double  *gdec;          /* Array of declinations (returned) */
     double  *gmag;          /* Array of magnitudes (returned) */
     int     *gtype;         /* Array of object types (returned) */
     int     nlog;           /* 1 for diagnostics */
{
  
  /*   int		name_given=0;	*/	/* Name or lists? */ 
  
  /* Request fields */
  char	survey[64], equinox[64], email[64];
  char	box[64], img[64], list[64];
  char	coords[1024], numbers[64];
  float sbox;
  int	i,j;
  char	*p;
  double mag;
  char  cs;
  /*Variables for getimage */
  int s;
  char buf[16384];
/*   char nbuf[16384]; */
  char nbuf2[16384];
  
  struct sockaddr_in sin;
  struct hostent *hp;
/*   struct servent *sp; */
  int val;
/*   //int i; */
  /*Variables for ptrans */
/*   float fnul; */
/*   char  sign; */
  int   decs;
  char	line[16384];
  char	block[16384];
  int   ns=0,llin=0;
  int 	lcnt = 0;
/*   int	scaling = 0; */
/*   char	lbuf[20]; */
/*   char	bigbuf[2880]; */
  int	nlset=0;
  char	c;
  int	fp=1;
  int   nl=0;
  float se,as;
  int ho,mi,dg,am;
  char  w1[1000],w2[1000],w3[1000];
  
  /* Initialize all fields to zero */
  survey[0]  = 0;
  equinox[0] = 0;
  email[0] = 0;
  box[0] = 0;
  img[0] = 0;
  list[0] = 0;
  coords[0] = 0;
  numbers[0] = 0;
  if(!(classd==-1 || classd==0 || classd==1 || classd==2)) {
    if(nlog) printf(" APM Class type %d not valid\n",classd);
    return(0);
  }

  /* make mag1 always the smallest magnitude */
  if (mag2 < mag1) {
    mag = mag2;
    mag2 = mag1;
    mag1 = mag;
  }
  

  /* Default output file */
  
/*   //out_file = "image.ps" ; */

  
  /* Loop over other arguments */	     
/*   //strcpy(email,"EMAIL"); */
/*   //strcat(email,"null"); */
  cra=cra/15.;
/*   //printf("cra %f\n",cra); */
  ho=(int)(cra);
  mi=(int)((cra-ho)*60);
  se=((cra-ho)*60-mi)*60;
/*   //printf(" ar %d %d %f \n",ho,mi,se); */
  dg=(int)(fabs(cdec));
  am=(int)((fabs(cdec)-dg)*60);
  as=((fabs(cdec)-dg)*60-am)*60;
/*   //printf(" dec %d %d %f\n",dg,am,as); */
  if(cdec<0) cs='-';
  else cs='+';
/*   //  strcpy(coords,"RADEC=3 42 10.0 +0 3 48"); */
  sprintf(coords,"RADEC=%2d %2d %4.1f %c%02d %2d %4.1f",ho,mi,se,cs,dg,am,as);
/*   //printf(" RADEC %s\n",coords); */
  strcpy(survey, "CAT=");
  strcat(survey, "poss1");
  
/*   //  strcpy(equinox, "EQUINOX="); */
/*   //  strcat(equinox, "j2000"); */
  sprintf(equinox,"EQUINOX=j2000");
  
  dra=dra*cos(cdec/180*3.1415);
  sbox= dra>ddec?dra:ddec;
  strcpy(box, "BOX=");
  strcat(box, "15");
  sprintf(box,"BOX=%f",sbox*60*2);
  
  strcpy(numbers, "NUMBERS=");
  strcat(numbers, "n");
  
  strcpy(list, "LIST=");
  strcat(list, "tmp_apm_request.list");
/*   out_file = "tmp_apm_request.list" ; */
/*   *out_file=0; */
  
/*   printf(" APMREAD OPTIONS : <<coor %s sur %s equin %s list %s box %s num %s img %s emial %s>>\n",coords,survey,equinox,list,box,numbers,img,email); */
  
  
  p= form_request(coords, survey, equinox, list, box, numbers,
		  img, email);
  
  strcpy(req_str, p);
  
  /*    fprintf(stderr, "%s\n", req_str); */
  /* Comienza getimage*/
/*   //------------------------------------------- */
  hp = gethostbyname(APMCAT);
  if (hp == NULL) {
    return ERROR;
  }
  /*
   *  Create an IP-family socket on which to make the connection
   */
  
  s = socket(hp->h_addrtype, SOCK_STREAM, 0);
  if (s < 0) {
    printf(" ERROR: Unable to create socket\n");
    return ERROR;
  }
    
  sin.sin_family = hp->h_addrtype;
  memcpy(&sin.sin_addr, hp->h_addr, hp->h_length);
  sin.sin_port = htons(PORT);
/*   //printf(" hp information: name %s alias %s add %d len %d lis %s\n",hp->h_name,hp->h_aliases[0],hp->h_addrtype,hp->h_length,hp->h_addr); */
 
  /*
   *  Connect to that address...
   */
  
  if ((val=connect(s, &sin, sizeof(sin))) < 0) {
/*     printf(" val %d\n",val); */
    printf(" ERROR: Unable to connect server APMCAT %s\n",APMCAT);
    return ERROR;
  }
  printf(" APMCAT server connected. Waiting for reply\n");
  
  /* Send first part of request */
  strcpy(buf, "POST /apmcatbin/post-query HTTP/1.0\012");
  
  if (send(s, buf, strlen(buf),0) < 0) {
     printf(" ERROR: While trying to send package to APMCAT\n");
     return ERROR;
  }  
  strcpy(buf, "User_Agent: SkyView Image Selector\012");
  strcat(buf, "Content-type: application/x-www-form-urlencoded\012");
  sprintf(nbuf2, "Content-length: %d\012\012", strlen(req_str));
  
  strcat(buf, nbuf2);
  strcat(buf, req_str);
/*   //printf(" BUFFER <<%s>>\n",buf); */
  
  /* Send request contents */
  if (send(s, buf, strlen(buf), 0) < 0) {
     printf(" ERROR: While trying to send package to APMCAT\n");
     return ERROR;
  }
  printf(" Requesting data... Please wait.\n");
  /* printf("%s",buf); */
  
  /* Now get back the data */
/*   //  if (ptrans(s) < 0)  return ERROR; */

/*   //Comienza ptrans */
/*   //-------------------------------------------- */
  
  /* Look for a double newline to signal the beginning of the data */
  while (1) {
    if (getbuf(s, &c, 1) < 0) 
      return ERROR;
    if (c == '\012') 
      {
	if (nlset)
	  break;
	else
	  nlset = 1;
      } else nlset = 0;
  }
  
  if (out_file) if (*out_file) 
    {
      if (info_req)
	fp = creat(out_file, 0644);
      else
	{
	  fp = creat(out_file,0644);
	}
      if (fp <=0) 
	{
	  printf("ERROR: Unable to open output file %s\n", out_file);
	  perror("  creat");
	}
    }
  
/*   //nl=0;ns=0;llin=0; */
  lcnt=getbuf(s, block, 1);
  printf(" lcnt %d BLOCK <<%s>>\n",lcnt,block); 
  while(1) {
    lcnt=getbuf(s, block, 2880);
    if(strstr(block,"Invalid parameters")!=NULL) {
      printf(" APMREAD: Invalid parameters\n");
/*       //return(ns); */
    }
    printf(" lcnt %d BLOCK <<%s>>\n",lcnt,block); 
    if (lcnt <= 0) break;
    for(i=0;i<lcnt;i++) { 
      printf(" bloque <<%c>> %d\n",block[i],i); 
      if(block[i]==24) return(ns);
      
      if(block[i]=='\n' || block[i]==EOF) { 
	nl++; 

 	if(nl>9) {  
	  printf(" LINEA %d <<%s>>\n",nl,line); 
 	  line[llin]='\n';  
	  LeeWord(line,1,w1);LeeWord(line,2,w2);LeeWord(line,3,w3);
 	  gra[ns]=(atof(w1)+atof(w2)/60.+atof(w3)/3600.)*15;  
	  LeeWord(line,4,w1);LeeWord(line,5,w2);LeeWord(line,6,w3);
 	  if(w1[0]=='-') decs=-1;  
 	  else decs=+1;  
/* 	  //	  printf(" w1 %s w2 %s w3 %s dec %d\n",w1,w2,w3,decs); */
	  gdec[ns]=decs*(fabs(atof(w1))+atof(w2)/60.+atof(w3)/3600.);  
	  LeeWord(line,7,w1);
	  gmag[ns]=atof(w1);
	  LeeWord(line,8,w1);
	  gtype[ns]=(int)atof(w1);
	  LeeWord(line,19,w1);
	  gnum[ns]=atof(w1);
	  if(nlog) printf(" star %d ar %f dec %f mag %f num %f type %d \n",
                   ns,gra[ns],gdec[ns],gmag[ns],gnum[ns],gtype[ns]);  
/* 	  //ns++; */
	  
 	  if(gtype[ns]==classd) {
	    if((gmag[ns]>mag1 && gmag[ns]<mag2) || mag1==mag2) ns++;  
	  }
	  if(ns>=nstarmax) return (ns);
 	}  
	llin=0; 
	for(j=0;j<16384;j++) line[j]=0; 
      } 
      else { 
/* 	//printf(" llin %d\n",llin);  */
	line[llin]=block[i]; 
/* 	//printf(" line <<%s>>\n",line);  */
	llin++; 
/* 	//printf(" Pasa aqui\n");  */
      } 
/*        //printf(" Siguiente\n");  */
    } 
/*     //printf(" Paso por aqui\n"); */
/*     //write(fp, block, lcnt); */

/*     //nl++; */
  }
  return ns;
  

/*   //getimage(); */
  
  
}

char *form_request(coords, survey, equinox, list, box, numbers,
		   img, email)
     
     char	*coords, *survey, *equinox, *list, *box, *numbers;
     char	*img, *email;
     
{
  static char buf[2048];
/*   char	lbuf[128]; */
/*   char	*t = buf; */
  
  
  /* Concatenate argument fields together.  Add in defaults
   * for arguments not specified.
   */
  strcpy(buf, mo_code(coords));
  
  if ( *list && *img ) {
    printf("\n\n Please specify either an image OR a list\n\n");
    exit(1);
  } 
  
  strcat(buf, "&");
  if (*survey)
    {
      
      strcat(buf, mo_code(survey));
    } else 
      {
	strcat(buf, "SURVEY=poss1");
      }
  
  strcat(buf, "&");
  if (*list)
    {
      strcat(buf, mo_code(list));
    } else 
      {
	strcat(buf, "LIST=off");
      }
  
  strcat(buf, "&");
  if (*equinox)
    {
      strcat(buf, mo_code(equinox));
    } else 
      {
	strcat(buf, "EQUINOX=b1950");
      }
  
  strcat(buf, "&");
  if (*box)
    {
      strcat(buf, mo_code(box));
    } else 
      {
	strcat(buf, "BOX=5");
      }
  
  strcat(buf, "&");
  if (*numbers)
    {
      strcat(buf, mo_code(numbers));
    } else 
      {
	strcat(buf, "NUMBERS=n");
      }
  
  strcat(buf, "&");
  if (*img)
    {
      strcat(buf, mo_code(img));
    } else 
      {
	strcat(buf, "PS=file.ps");
      }
  
  strcat(buf, "&");
  if (*email)
    {
      strcat(buf, mo_code(email));
    } else 
      {
	strcat(buf, "EMAIL=null");
      }
  
  strcat(buf, "\013\012\012");
  printf(" buf <<%s>>\n",buf);
  
  return buf;
}


/* Send request and receive data.
 */
/* int getimage() */
     
/* { */
/*   int s, n; */
/*   char buf[16384]; */
/*   char nbuf[16384]; */
/*   char nbuf2[16384]; */
  
/*   struct sockaddr_in sin; */
/*   struct hostent *hp; */
/*   struct servent *sp; */
/*   int val; */
/*   int i; */
  
/*   hp = gethostbyname(APMCAT); */
/*   if (hp == NULL) { */
/*     return ERROR; */
/*   } */
/*    */
/*    *  Create an IP-family socket on which to make the connection */
/*    */ 
  
/*   s = socket(hp->h_addrtype, SOCK_STREAM, 0); */
/*   if (s < 0) { */
/*     return ERROR; */
/*   } */
  
  
/*   sin.sin_family = hp->h_addrtype; */
/*   memcpy(&sin.sin_addr, hp->h_addr, hp->h_length); */
/*   sin.sin_port = htons(PORT); */
  
/*    */
/*    *  Connect to that address... */
/*    */ 
  
/*   if ((val=connect(s, &sin, sizeof(sin))) < 0) { */
/*     return ERROR; */
/*   } */
  
/* Send first part of request */ 
/*   strcpy(buf, "POST /apmcatbin/post-query HTTP/1.0\012"); */
  
/*   if (send(s, buf, strlen(buf),0) < 0) */
/*     return ERROR; */
  
/*   strcpy(buf, "User_Agent: SkyView Image Selector\012"); */
/*   strcat(buf, "Content-type: application/x-www-form-urlencoded\012"); */
/*   sprintf(nbuf2, "Content-length: %d\012\012", strlen(req_str)); */
  
/*   strcat(buf, nbuf2); */
/*   strcat(buf, req_str); */
  
/* Send request contents */ 
/*   if (send(s, buf, strlen(buf), 0) < 0) */
/*     return ERROR; */
/* printf("%s",buf); */ 
  
/* Now get back the data */ 
/*   if (ptrans(s) < 0) */
/*     return ERROR; */
/* } */

char	*mo_code(string)
     
     char	*string;
     
{
  /* Perform Mosaic encoding.  Only alphanumerics are
   * unchanged.  Spaces are replaced by +.  All others by
   * %xx where xx is the Hex code for the value.
   * Note this routine uses a static pointer so the value
   * should be copied before this routine is called again.
   */
  static char buf[2048];
  char	*s,*t;
  char	c;
  
  s = string;
  t = buf;
  do {
    
    c = *s;
    *t++ = *s++;
  }
  
  while (c != '=');
  
  c = *s;
  while (c) 
    {
      if (c == ' ') *t++ = '+';
      
      /* Assume ASCII sequencing */
      else if ( !  ((c >= 'a' && c <= 'z') ||
		    (c >= 'A' && c <= 'Z') ||
		    (c >= '0' && c <= '9')
		    )) {
	sprintf(t, "%%%2x", c);
	t += 3;
      } else 
	{
	  *t = c;
	  t++;
	}
      
      
      c = *++s;
    }
  
  *t = 0;
  return buf;
}

/* Read from input until we can fill the buffer */
int getbuf(s, rbuf, qlen)
     
     int s;
     char	*rbuf;
     int qlen;
{
  
  
  static int left=0;
  static char bigbuf[32768];
  int	n;
  int	len = qlen;
  int     i;
  
  while (len > left) {
    n = recv(s, bigbuf+left, sizeof(bigbuf)-left, 0);
    if (n <= 0) {
      if (left <= 0)
	return( -1);
      else {
	len = left;
	break;
      }
    }
    else 
      {
	left = left + n;
      }
  }
  memcpy(rbuf, bigbuf, len);
  left = left - len;
  if (left > 0)
    for (i = 0; i < left; i++)
      bigbuf[i] = bigbuf[i+len];
  /*	        memmove(bigbuf, bigbuf+len, left);     */
  return len;
}

/* Transfer data from socket to standard output
 */
/* int ptrans(s) */
     
/*      int	s; */
     
/* { */
/*   char	line[16384]; */
/*   int 	lcnt = 0; */
/*   int	scaling = 0; */
/*   char	lbuf[20]; */
/*   char	bigbuf[2880]; */
/*   int	nlset=0; */
/*   char	c; */
/*   int	fp=1; */
/*   int   nl; */
  
/* Look for a double newline to signal the beginning of the data */ 
/*   while (1) { */
/*     if (getbuf(s, &c, 1) < 0)  */
/*       return ERROR; */
/*     if (c == '\012')  */
/*       { */
/* 	if (nlset) */
/* 	  break; */
/* 	else */
/* 	  nlset = 1; */
/*       } else nlset = 0; */
/*   } */
  
/*   if (out_file) if (*out_file)  */
/*     { */
/*       if (info_req) */
/* 	fp = creat(out_file, 0644); */
/*       else */
/* 	{ */
/* 	  fp = creat(out_file,0644); */
/* 	} */
/*       if (fp <=0)  */
/* 	{ */
/* 	  printf("ERROR: Unable to open output file %s\n", out_file); */
/* 	  perror("  creat"); */
/* 	} */
/*     } */
  
/*   nl=0; */
/*   lcnt=getbuf(s, line, 1); */
/*   while(1) { */
/*     lcnt=getbuf(s, line, 2880); */
/*     if (lcnt <= 0) break; */
/*     write(fp, line, lcnt); */
/*   } */
/*   return 0; */
/* } */

