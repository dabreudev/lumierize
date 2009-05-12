#define DEBUG 0
/*  This file, cfileio.c, contains the low-level file access routines.     */

/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.  Users shall not, without prior written   */
/*  permission of the U.S. Government,  establish a claim to statutory     */
/*  copyright.  The Government and others acting on its behalf, shall have */
/*  a royalty-free, non-exclusive, irrevocable,  worldwide license for     */
/*  Government purposes to publish, distribute, translate, copy, exhibit,  */
/*  and perform such material.                                             */
#include "modulos.h"

/*  #include <string.h> */
/*  #include <stdlib.h> */
/*  #include <math.h> */
#include <ctype.h>
/*  #include <stddef.h>  */
/*  #include "fitsio2.h" */
/*  #include "group.h" */

/*Aqui empieza un trozo de fitsio2.h */
#ifndef _FITSIO2_H
#define _FITSIO2_H

#define DBUFFSIZE 28800 /* size of data buffer in bytes */
 
#define NIOBUF  25       /* number of IO buffers to create */
#define IOBUFLEN 2880    /* size in bytes of each IO buffer */
#define MINDIRECT 8640   /* minimum size for direct reads and writes */
                         /* MINDIRECT must have a value >= 8640 */
 
#define NATIVE             0 /* a generic machine that uses IEEE formats */
#define ULTRIX             1
#define ALPHA_OSF          2
#define VAXVMS             3
#define ALPHAVMS           4
#define IBMPC              5
#define CRAY               6
 
#define GFLOAT             1
#define IEEEFLOAT          2
#define DATA_UNDEFINED -1      
#define REPORT_EOF 0                   
#define maxvalue(A,B) ((A) > (B) ? (A) : (B))       
#endif                
/* Hasta aqui */
                                                                                                       


#define MAX_PREFIX_LEN 20  /* max length of file type prefix (e.g. 'http://') */
#define MAX_DRIVERS 20     /* max number of file I/O drivers */

typedef struct    /* structure containing pointers to I/O driver functions */ 
{   char prefix[MAX_PREFIX_LEN];
    int (*init)(void);
    int (*shutdown)(void);
    int (*setoptions)(int option);
    int (*getoptions)(int *options);
    int (*getversion)(int *version);
    int (*checkfile)(char *urltype, char *infile, char *outfile);
    int (*open)(char *filename, int rwmode, int *driverhandle);
    int (*create)(char *filename, int *drivehandle);
    int (*truncate)(int drivehandle, long size);
    int (*close)(int drivehandle);
    int (*remove)(char *filename);
    int (*size)(int drivehandle, long *size);
    int (*flush)(int drivehandle);
    int (*seek)(int drivehandle, long offset);
    int (*read)(int drivehandle, void *buffer, long nbytes);
    int (*write)(int drivehandle, void *buffer, long nbytes);
} fitsdriver;

fitsdriver driverTable[MAX_DRIVERS];  /* allocate driver tables */

/*  int need_to_initialize = 1; */    /* true if CFITSIO has not been initialized */
/*  int no_of_drivers = 0;        */  /* number of currently defined I/O drivers */


/*--------------------------------------------------------------------------*/
int ffopen2(fitsfile **fptr,      /* O - FITS file pointer                   */ 
           const char *name,     /* I - full name of file to open           */
           int mode,             /* I - 0 = open readonly; 1 = read/write   */
           int *status)          /* IO - error status                       */
/*
  Open an existing FITS file with either readonly or read/write access.
*/
{
    int  driver, hdutyp, slen, writecopy;
    long filesize;
    int  handle;
    char urltype[MAX_PREFIX_LEN], infile[FLEN_FILENAME], outfile[FLEN_FILENAME];
    char origurltype[MAX_PREFIX_LEN], extspec[FLEN_FILENAME];
    char rowfilter[FLEN_FILENAME];
    char imagecolname[FLEN_VALUE], rowexpress[FLEN_FILENAME];
    char binspec[FLEN_FILENAME], colspec[FLEN_FILENAME];
    char histfilename[FLEN_FILENAME];
    char filtfilename[FLEN_FILENAME];

    char *url;

    int found_end=0;
    int nextkey=0;
    char  keyname[FLEN_KEYWORD], value[FLEN_VALUE], comm[FLEN_COMMENT];
    int nspace=0;

    if(DEBUG) printf("Entro en ffopen2\n");


    if (*status > 0)
        return(*status);

    *fptr = 0;              /* initialize null file pointer */
    writecopy = 0;  /* have we made a write-able copy of the input file? */

/*      if (need_to_initialize)    */        /* this is called only once */
    /* Supongo que si que hay que inicializar */
    if(DEBUG) printf("Aqui inicializo \n");
    *status = fits_init_cfitsio();

    if(DEBUG) printf("Acabo de inicializar\n");

    if (*status > 0)
        return(*status);

    url = (char *) name;
    while (*url == ' ')  /* ignore leading spaces in the filename */
        url++;

    if (*url == '\0')
    {
        ffpmsg("Name of file to open is blank. (ffopen)");
        return(*status = FILE_NOT_OPENED);
    }

        /* parse the input file specification */
    ffiurl(url, urltype, infile, outfile, extspec,
              rowfilter, binspec, colspec, status);

    if(DEBUG) printf(" infile %s, urltype %s, outfile %s, extspec %s, rowfilter %s, binspec %s, colspec %s\n",infile,urltype,outfile,extspec,rowfilter,binspec,colspec);

    
    if(DEBUG) printf("Ya se que tipo es\n");

    if (*status > 0)
    {
        ffpmsg("could not parse the input filename: (ffopen)");
        ffpmsg(url);
        return(*status);
    }

    imagecolname[0] = '\0';
    rowexpress[0] = '\0';

    histfilename[0] = '\0';
    filtfilename[0] = '\0';
    

    if(DEBUG) printf("Aqui ya ha hecho cositas\n");

    /*-------------------------------------------------------------------*/
    /* check if this same file is already open, and if so, attach to it  */
    /*-------------------------------------------------------------------*/

    *status = urltype2driver(urltype, &driver);
    if (*status > 0)
    {
        ffpmsg("could not find driver for this file: (ffopen)");
        ffpmsg(urltype);
        ffpmsg(url);
        return(*status);
    }

    /*-------------------------------------------------------------------
        deal with all those messy special cases which may require that
        a different driver be used:
            - is disk file compressed?
            - are ftp: or http: files compressed?
            - has user requested that a local copy be made of
              the ftp or http file?
      -------------------------------------------------------------------*/

    if (driverTable[driver].checkfile)
    {
      if(DEBUG) printf("Se ha metido en un lio\n");
      
      strcpy(origurltype,urltype);  /* Save the urltype */
      
      /* 'checkfile' may modify the urltype, infile and outfile strings */
      *status =  (*driverTable[driver].checkfile)(urltype, infile, outfile);

        if (*status)
        {
            ffpmsg("checkfile failed for this file: (ffopen)");
            ffpmsg(url);
            return(*status);
        }

        if (strcmp(origurltype, urltype))  /* did driver changed on us? */
        {
            *status = urltype2driver(urltype, &driver);
            if (*status > 0)
            {
                ffpmsg("could not change driver for this file: (ffopen)");
                ffpmsg(url);
                ffpmsg(urltype);
                return(*status);
            }
        }
    }
    
    if(DEBUG) printf("Antes de otro lio\n");

    /* call appropriate driver to open the file */
    if (driverTable[driver].open)
    {
      if(DEBUG) printf("Se ha metido en otro\n"); 
      
      *status =  (*driverTable[driver].open)(infile, mode, &handle);
        if (*status > 0)
        {
            ffpmsg("failed to find or open the following file: (ffopen)");
            ffpmsg(url);
            return(*status);
       }
    }
    else
    {
        ffpmsg("cannot open an existing file of this type: (ffopen)");
        ffpmsg(url);
        return(*status = FILE_NOT_OPENED);
    }

    if(DEBUG) printf("Esto es antes de ver el tamanio\n");

        /* get initial file size */
    *status = (*driverTable[driver].size)(handle, &filesize);
    if (*status > 0)
    {
        (*driverTable[driver].close)(handle);  /* close the file */
        ffpmsg("failed get the size of the following file: (ffopen)");
        ffpmsg(url);
        return(*status);
    }
    if(DEBUG) printf("Essto es antes de alocatear\n");

        /* allocate fitsfile structure and initialize = 0 */
    *fptr = (fitsfile *) calloc(1, sizeof(fitsfile));

    if(DEBUG) printf("Essto es despues de alocatear\n");    

    if (!(*fptr))
    {

      if(DEBUG) printf("Esto es cerrar el handle\n");

        (*driverTable[driver].close)(handle);  /* close the file */
        ffpmsg("failed to allocate structure for following file: (ffopen)");
        ffpmsg(url);
        return(*status = MEMORY_ALLOCATION);
    }

    if(DEBUG) printf("Alocateo estrucutra\n");

        /* allocate FITSfile structure and initialize = 0 */
    (*fptr)->Fptr = (FITSfile *) calloc(1, sizeof(FITSfile));



    if (!((*fptr)->Fptr))
    {
      
      if(DEBUG) printf("Vuelvo a cerrar!\n");

        (*driverTable[driver].close)(handle);  /* close the file */
        ffpmsg("failed to allocate structure for following file: (ffopen)");
        ffpmsg(url);
        free(*fptr);
        *fptr = 0;       
        return(*status = MEMORY_ALLOCATION);
    }

    if(DEBUG) printf("Calcula longitued \n");

    slen = strlen(url) + 1;
    slen = maxvalue(slen, 32); /* reserve at least 32 chars */ 
    ((*fptr)->Fptr)->filename = (char *) malloc(slen); /* mem for file name */

    if ( !(((*fptr)->Fptr)->filename) )
    {

      

        (*driverTable[driver].close)(handle);  /* close the file */
        ffpmsg("failed to allocate memory for filename: (ffopen)");
        ffpmsg(url);
        free((*fptr)->Fptr);
        free(*fptr);
        *fptr = 0;              /* return null file pointer */
        return(*status = MEMORY_ALLOCATION);
    }

    if(DEBUG) printf("No se que de tablas\n");


        /* store the parameters describing the file */
    ((*fptr)->Fptr)->filehandle = handle;        /* file handle */
    ((*fptr)->Fptr)->driver = driver;            /* driver number */
    strcpy(((*fptr)->Fptr)->filename, url);      /* full input filename */
    ((*fptr)->Fptr)->filesize = filesize;        /* physical file size */
    ((*fptr)->Fptr)->logfilesize = filesize;     /* logical file size */
    ((*fptr)->Fptr)->writemode = mode;           /* read-write mode    */
    ((*fptr)->Fptr)->datastart = DATA_UNDEFINED; /* unknown start of data */
    ((*fptr)->Fptr)->curbuf = -1;            /* undefined current IO buffer */
    ((*fptr)->Fptr)->open_count = 1;      /* structure is currently used once */
    ((*fptr)->Fptr)->validcode = VALIDSTRUC; /* flag denoting valid structure */

    ffldrc(*fptr, 0, REPORT_EOF, status);     /* load first record */

    if(DEBUG) printf("Aqui acaba lo esstablecido\n");

/*      if (ffrhdu(*fptr, &hdutyp, status) > 0)  */ /* determine HDU structure */
/*      { */
/*          ffpmsg( */
/*            "ffopen could not interpret primary array header of file: "); */
/*          ffpmsg(url); */

/*          if (*status == UNKNOWN_REC) */
/*             ffpmsg("This does not look like a FITS file."); */

/*          ffclos(*fptr, status); */
/*          *fptr = 0;    */           /* return null file pointer */
/*          return(*status); */
/*      } */

    /* Aqui meto yo lo nuevo */
 

    
/*      ffpinit(fptr, status);       */      /* initialize the primary array */
    
    hdutyp = 0;
 
    ((*fptr)->Fptr)->hdutype = IMAGE_HDU; /* primary array or IMAGE extension  */
    ((*fptr)->Fptr)->headend = ((*fptr)->Fptr)->logfilesize;  /* set max size */
     
    if(DEBUG) printf("HAce chirivitas\n");
    
    for (; !found_end; nextkey++)      {
      /* get next keyword */
      if(DEBUG) printf(" 666 Me voy por %d\n",nextkey); 
      if (ffgkyn(*fptr, nextkey, keyname, value, comm, status) > 0)
	{
	  fits_report_error(stdout,*status);
	  
	  if(DEBUG)   printf("ERROR Me voy por %d\n",nextkey); 
	  if (*status == KEY_OUT_BOUNDS)
	    {
	      found_end = 1;  /* simply hit the end of the header */
	    }
	  else          
	    {
	      ffpmsg("Failed to find the END keyword in header (ffgphd).");
	    }
	}
      else   /* got the next keyword without error */
	{
	  
	  if(DEBUG) printf("antes str\n");
	  if(DEBUG) printf("Teiene que hacer esto keyname <%s> value  <%s> com  <%s>\n",keyname,value,comm);
	  
	  if (!strcmp(keyname, "END"))
	    found_end = 1;
	  
	  else if (!keyname[0] && !value[0] && !comm[0])
	    nspace = nspace + 1;  /* this is a blank card in the header */
	  
	  else
	    nspace = 0;  /* reset count of blank keywords immediately
			    before the END keyword to zero   */
	  
	  if(DEBUG) printf("Ya no else\n");
	  
	}
    } 
    
    if(DEBUG) printf("Lupecito salio\n");

    if(DEBUG) printf("nspace %d\n",nspace);


    ((*fptr)->Fptr)->headend = ((*fptr)->Fptr)->nextkey  - (80 * (nspace + 1));

    if(DEBUG) printf("primero\n");


    ((*fptr)->Fptr)->datastart = ( (((*fptr)->Fptr)->nextkey - 80) / 2880 + 1) * 2880;
    ((*fptr)->Fptr)->headstart[ ((*fptr)->Fptr)->curhdu + 1] =
      ((*fptr)->Fptr)->datastart + 
      ( 0 + 2879) / 2880 * 2880;

    if(DEBUG) printf("algo heads\n");


    ((*fptr)->Fptr)->heapstart = 0;
    ((*fptr)->Fptr)->heapsize = 0;
    ((*fptr)->Fptr)->rowlength = 0;    /* rows have zero length */
    ((*fptr)->Fptr)->tfield = 0;       /* table has no fields   */
    ((*fptr)->Fptr)->tableptr = 0;     /* set a null table structure pointer */
    
    ((*fptr)->Fptr)->nextkey = ((*fptr)->Fptr)->headstart[ ((*fptr)->Fptr)->curhdu ];
    
    /*  compare the starting position of the next HDU (if any) with the size */
    /*  of the whole file to see if this is the last HDU in the file */

    ((*fptr)->Fptr)->lasthdu = 0;  /* no, not the last HDU */
    ((*fptr)->Fptr)->lasthdu = 1;  /* yes, this is the last HDU */
     
    
    if(DEBUG) printf("Y aqui acabo\n");

    
    /* ------------------------------------------------------------- */
    /* At this point, the input file has been opened. If outfile was */
    /* specified, then we have opened a copy of the file, not the    */
    /* original file so it is safe to modify it if necessary         */
    /* ------------------------------------------------------------- */

    if (*outfile)
        writecopy = 1;  



    return(*status);
}







