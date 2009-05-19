#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# TITLE: Script para recolectar simulaciones con LumFunc_GOYA
# CREATED: 20080325 dabreu
# Este script recolecta los resultados de una simulación y los pone en el
# fichero correspondiente.
#
# MODIFIED: 20080402 dabreu
# Añadidas dos columnas nuevas a los ficheros .lf
#
# Uso: LF_harvester.py -m method simulationParameterFile
#

"""Este script recolecta el resultado de una simulación"""

import sys
from optparse import OptionParser
import os
import catg

simulationsDir=os.getenv("LF_SIM")

#Función para leer los parámetros de la simulación
def readParams(file):
    """readParams(file) -> devuelve un diccionario con los datos del fichero de
                        parámetros en un diccionario y los valores como
			strings"""
    parametros={}
    f=open(file,'r')
    lineas=f.readlines()
    f.close()
    for linea in lineas:
        if linea[0]=='#':
            continue
        else:
            trozos=linea.split('#')[0]
	    trozos=trozos.split()
	    if len(trozos[1:])>1: #parámetro con más de un valor
	        parametros[trozos[0]]=trozos[1:] #se devuelve como una lista
	    else:
	        parametros[trozos[0]]=trozos[1]

    return parametros

def sumarCatalogos(cat1,cat2):
    """sumarCatalogos(cat1,cat2) -> Devuelve un catálogo resultado de añadir
   				cat2 a cat1"""
    for i in cat1['llaves']:
        cat1[i]=list(cat1[i])
        for j in cat2[i]:
            cat1[i].append(j)
        try:
            cat1[i]=numpy.asarray(cat1[i])
        except:
            pass
	 
    return cat1

def main():
    usage = "usage: %prog [options] simulationParameterFile \nDefault values: out=stdout, verbose=1 method=VVmax"
    parser = OptionParser(usage)
    parser.add_option("-o", "--out", dest="outfilename", default="stdout",
                      help="Output to FILE [default: %default]", metavar="FILE")
    parser.add_option("-v", "--verbose", dest="verbose", default=1,
                      help="Set verbose level to INT [default: %default", metavar="INT")
    parser.add_option("-m", "--method", dest="type", default='VVmax',
                choices=["VVmax","STY","STY_wC","STY_errMag","STY_wC_errColor"],
                help="Simulation TYPE [default: %default]", metavar="TYPE")
    parser.set_default("outfilename","stdout")

    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Error in args. Type --help for help\n"
        sys.exit(1)

    simulationParameterFile=args[0]
    verbose=options.verbose
    method=options.type
      
    if options.outfilename=="stdout":
        output=sys.stdout
    else:
        output=open(options.outfilename, "w")
      
    parameterFile=simulationsDir+"/"+simulationParameterFile
    paramValues=readParams(parameterFile)

    niter=int(paramValues['NITER'])   
    tipo=paramValues['TYPE']
    H0=float(paramValues['H0'])
    OM=float(paramValues['OM'])
    OL=float(paramValues['OL'])
    Mstar=float(paramValues['MSTAR'])
    Phistar=float(paramValues['PHISTAR'])
    Alfa=float(paramValues['ALFA'])
    zlow=float(paramValues['ZLOW'])
    zup=float(paramValues['ZUP'])
    zerr=float(paramValues['ZERR'])
    zderr=float(paramValues['ZDERR'])
    mlim=float(paramValues['MLIM'])
    deltamag=float(paramValues['DELTAMAG'])
    merr=float(paramValues['MERR'])
    mderr=float(paramValues['MDERR'])
    color=float(paramValues['C'])
    colord=float(paramValues['CD'])
    cerr=float(paramValues['CERR'])
    cderr=float(paramValues['CDERR'])
    density=float(paramValues['DENSITY'])
    nobjects=[float(i) for i in paramValues['NOBJECTS']]

    for numObjects in nobjects:
        if verbose: output.write('Dentro del bucle. '+str(int(numObjects))+'\n')
        area=numObjects/density

        #Names of files
        simulationName="_H0"+str(H0)+"_OM"+str(OM)+"_OL"+str(OL)
        simulationName+="_mstar"+str(Mstar)+"_phistar"+str(Phistar)
        simulationName+="_alfa"+str(Alfa)
        simulationName+="_zlow"+str(zlow)+"_zup"+str(zup)
        simulationName+="_zerr"+str(zerr)+"_zderr"+str(zderr)
        simulationName+="_mlim"+str(mlim)+"_deltm"+str(deltamag)
        simulationName+="_merr"+str(merr)+"_mderr"+str(mderr)
        simulationName+="_c"+str(color)+"_cd"+str(colord)
        simulationName+="_cerr"+str(cerr)+"_cderr"+str(cderr)
        simulationName+="_area"+str(area)+"/"


        resultsDir=simulationsDir+"/Results/"+method+"/cosecha"+simulationName+"/"

        if not(os.access(resultsDir,os.F_OK)):
            os.makedirs(resultsDir)

        finalName=resultsDir+method+"_numObjects_"+str(int(numObjects))+".cat"

        allData=[]
        for iter in range(niter):

            computedFile=simulationsDir+"/ComputedLF/"+method+"/simlf_m"+simulationName+"/"+str(iter)+".lf"

            # done files
            doneFileName=simulationsDir+"/LockFiles/"+method+"/simlf_m"+simulationName+"/"+str(iter)+".done"
            if os.access(doneFileName,os.F_OK): 
                toReadFile=open(computedFile,'r')
                lines=toReadFile.readlines()
                toReadFile.close()
                try:
                    data=catg.rcat(computedFile)
                    allData.append(data)
                except:
                    output.write("Big problems with "+computedFile+"\n")

            else:
                output.write("No se encontró a "+ doneFileName+"\n") 

        # All the files are now in allData
        finalCat=allData[0]
        allData.remove(allData[0])
        for data in allData:
            finalCat=sumarCatalogos(finalCat,data)

        catg.writecat(finalCat, finalName, overwrite=1)
        output.write("Final cat in "+finalName+"\n")

    if verbose: output.write("The harvest is ready.\n")   

    sys.exit()
      
if __name__ == "__main__":
    main()
