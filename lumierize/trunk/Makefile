OPERA_DIR = /home/cgarcia/OPERA
BIN_OPERA_DIR = /home/cgarcia/bin
LIB_OPERA_DIR = /home/cgarcia/lib


SHELL = /bin/sh
INSTALL = install
ARFLAGS = cr
RANLIB = $(AR) -s
INC_OPERA_DIR = $(OPERA_DIR)/include
PRG_OPERA_DIR = $(OPERA_DIR)/C_SOURCE
MOD_OPERA_DIR = $(OPERA_DIR)/modulos
#CFLAGS = -I$(GSL_DIR_INC) -I$(PGPLOT_DIR) -I$(CBUTTON_DIR) -I$(FITSIO_DIR) -I$(WCSTOOLS_DIR) -I$(INC_OPERA_DIR) -Wall -O2 -g -pg
CFLAGS = -I$(GSL_DIR_INC) -I$(PGPLOT_DIR) -I$(CBUTTON_DIR) -I$(FITSIO_DIR) -I$(WCSTOOLS_DIR) -I$(INC_OPERA_DIR) -Wall -O3 -g -march=pentiumpro 


export CFLAGS
export OPERA_DIR
export BIN_OPERA_DIR
export LIB_OPERA_DIR


redhat:	LDFLAGS = -L$(GSL_DIR_LIB)  -L$(PGPLOT_DIR)  -L$(FITSIO_DIR) -L$(CBUTTON_DIR) -L/usr/X11R6/lib/ -L$(WCSTOOLS_DIR) -L$(LIB_OPERA_DIR)  -lopera -lcbutton -lcpgplot -lpgplot -lg2c -lX11 -lcfitsio  -lm  -lwcs -lgsl -lgslcblas
	export LDFLAGS
redhat:   all 

linux: LDFLAGS = -L$(GSL_DIR_LIB)  -L$(PGPLOT_DIR)  -L$(FITSIO_DIR) -L$(CBUTTON_DIR) -L/usr/X11/lib/ -L$(WCSTOOLS_DIR) -L$(LIB_OPERA_DIR)  -lopera -lcbutton -lcpgplot -lpgplot -lf2c -lX11 -lcfitsio  -lm  -lwcs -lgsl -lgslcblas
	export linux
linux: all

digital: LDFLAGS = -L$(GSL_DIR_LIB)  -L$(PGPLOT_DIR)  -L$(FITSIO_DIR) -L$(CBUTTON_DIR) -L$(WCSTOOLS_DIR) -L$(LIB_OPERA_DIR)  -lopera -lcbutton -lcpgplot -lpgplot -lcfitsio  -lX11 -lm  -lfor -lUfor -lwcs -lgsl -lgslcblas
	export LDFLAGS
digital: all

sun: LDFLAGS = -L$(GSL_DIR_LIB)  -L$(PGPLOT_DIR)  -L$(FITSIO_DIR) -L$(CBUTTON_DIR) -L$(WCSTOOLS_DIR) -L/usr/X11R6/lib/ -L$(LIB_OPERA_DIR)  -lopera -lcbutton -lcpgplot -lpgplot -lcfitsio  -lX11 -lm  -lM77 -lF77 -lf77compat -lsunmath_mt  -lwcs -lgsl -lgslcblas
	export LDFLAGS
sun: all

solaris: LDFLAGS = -L$(GSL_DIR_LIB)  -L$(PGPLOT_DIR)  -L$(FITSIO_DIR) -L$(CBUTTON_DIR) -L$(WCSTOOLS_DIR) -L/usr/X11R6/lib/ -L$(LIB_OPERA_DIR)  -lopera -lcbutton -lcpgplot -lpgplot -lcfitsio  -lX11  -lM77 -lF77  -lsunmath_mt  -lm -lwcs -lgsl -lgslcblas
	export LDFLAGS
solaris: all

sparc: LDFLAGS = -L$(GSL_DIR_LIB)  -L$(PGPLOT_DIR)  -L$(FITSIO_DIR) -L$(CBUTTON_DIR) -L$(WCSTOOLS_DIR) -L$(LIB_OPERA_DIR)  -lopera -lcbutton -lcpgplot -lpgplot -lcfitsio  -lX11 -lm  -lfor -lUfor -lwcs -lgsl -lgslcblas
	export LDFLAGS
sparc: all

all: 
	@-$(MAKE) -s -C $(MOD_OPERA_DIR)
	@-$(MAKE) -s -C $(PRG_OPERA_DIR)
install: all 
	@-$(MAKE)  install -s -C $(MOD_OPERA_DIR)
	@-$(MAKE)  install -s -C $(PRG_OPERA_DIR)

.PHONY: clean
clean:
	@-$(MAKE)  clean -s -C $(MOD_OPERA_DIR)
	@-$(MAKE)  clean -s -C $(PRG_OPERA_DIR)
