#ifndef FERMISEL_H
#define FERMISEL_H

#ifdef __cplusplus
extern "C" {
#endif


	/* Para trabajar con funciones de seleccion */
	struct fermifsel_M  {
		double  magcut;
		double  deltamag;
	};
	struct fermifsel_L  {
		double  fluxcut;
		double  deltaflux;
	};
#ifdef __cplusplus
}
#endif


#endif
