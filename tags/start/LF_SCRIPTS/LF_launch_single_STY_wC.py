#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# TITLE: Script para generar calcular la LF con STY_wC
# CREATED: 20080325 dabreu
# Este script lanza el cálculo de la LF con STY_wC
# 
# MODIFIED: 20080423 dabreu
# Lockfiles sctructure modified
#
# Uso: LF_launch_single_STY_wC.py -p simulationParameterFile numberOfIter
#

"""Script para calcualr la LF con STY_wC"""

import sys
from optparse import OptionParser
import os
import signal
import time
#import commands
import subprocess
import datetime

estatus=0
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

def main():

    #Para poder capturar el sigterm y borrar el lockfile
    def signal_handler(signal, frame):
        """signal_handler(signal, frame) -> Para poder capturar el sigterm y borrar
                                       el lockfile

           Para que se produzca la captura, la función debe ser llamada con:
           signal.signal(signal.SIGTERM, signal_handler)"""

        output.write("Recibí un SIGTERM (habló launch_single_STY_wC)\n")
        try:
            os.remove(lockFileName)
            output.write("El lockfile se borró correctamente\n")
        except:
            output.write("La petó al intentar borrar el lockfile\n")
            pass
        if PID:
            try:
                os.kill(PID,15)
                time.sleep(0.05)
                output.write(str(PID)+" killed\n")
            except:
                output.write("Problems killing "+str(PID)+"\n")
        else:
            output.write("No such process "+str(PID)+"\n")

        sys.exit(-15)

    signal.signal(signal.SIGTERM, signal_handler)

    def secureCompute(nombrePaso,doneFileName,lockFileName,comando,outFileName):
        """secureCompute(nombrePaso,doneFileName,lockFileName,comando,outFileName)
      
   			   nombrePaso = catálogo, VVmax o STY (solo para verbose)
   			   doneFileName = nombre del fichero de done
   			   lockFileName = nombre del fichero de lock
   			   comando = lo que se enviará si todo está correcto
   			   outFileName = nombre del fichero de salida de comando"""
   
        doneFile=0; lockFile=0
   
        doneFile=os.access(doneFileName,os.F_OK)
        if not(doneFile): # El paso no está hecho
   	    if verbose: output.write("No está hecho el "+nombrePaso+"\n")
   	    lockFile=os.access(lockFileName,os.F_OK)
   	    if not(lockFile): # No se está haciendo
   	        os.system("touch "+lockFileName)
	        lockFile=os.access(lockFileName,os.F_OK)
	        if verbose: output.write("lockFile: "+str(lockFile)+"\n")
   	        if verbose: output.write(lockFileName+"\n")            
   	        if verbose: output.write("No se está haciendo el "+nombrePaso+"\n")
   	        #(estatus, salida)=commands.getstatusoutput(comando)
   	        #if verbose: output.write(salida)
                p=subprocess.Popen(comando,shell=True)
                PID=p.pid
                if verbose: output.write("PID: "+str(PID)+"\n")
                estatus=p.wait()
                if verbose: output.write("Terminó.\n")
   	        if estatus == 0: # El paso terminó bien
   	            if verbose: output.write("\nEl estatus es 0\n")
   	            if(os.access(outFileName,os.F_OK)):
   		        os.system("touch "+doneFileName)
   		        os.remove(lockFileName)
   	            else: # aparentemente terminó bien, pero no hay outFile
   		        if verbose: output.write("No aparece el "+outFileName+"\n")
   		        os.remove(lockFileName)
   		        return 1
   	        else: # El paso termino mal
   	            if verbose: output.write("La petó el "+nombrePaso+"\n")
	            output.write("=== "+lockFileName+"\n")
   	            os.remove(lockFileName)
   	            return 1
   	    
            else: # El paso se está haciendo
   	        if verbose: output.write("El "+nombrePaso+" se está haciendo\n")
   	        return 1 # pasamos al siguiente paso del for
   	    
        else: # El paso está hecho
   	    if verbose: output.write("El "+nombrePaso+" ya está hecho\n")
   
        return 0 # Todo termió bien

    # Aquí empieza main
    usage = "usage: %prog [options] Iter \nDefault values: out=stdout, par=default.par, verbose=1"
    parser = OptionParser(usage)
    parser.add_option("-o", "--out", dest="outfilename",
                     help="Output to FILE", metavar="FILE")
    parser.add_option("-p", "--par", dest="param",
                     help="Parameter filename (from WORK dir)", metavar="FILE")
    parser.add_option("-v", "--verbose", dest="verbose",
                     help="Set verbose level to INT", metavar="INT")
    parser.set_default("outfilename","stdout")
    parser.set_default("param","default.par")
    parser.set_default("verbose",1)
   
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Error in args. Type --help for help\n"
        sys.exit(1)

    simulationNumber=args[0]
    simulationParameterFile=options.param
    verbose=options.verbose
      
    if options.outfilename=="stdout":
        output=sys.stdout
    else:
        output=open(options.outfilename, "w")
      
    parameterFile=simulationsDir+"/"+simulationParameterFile
    paramValues=readParams(parameterFile)

    niter=int(paramValues['NITER'])   
    tipo=paramValues['TYPE'] # "m" or "l", but only "m" works
    Mstar=float(paramValues['MSTAR'])
    Phistar=float(paramValues['PHISTAR'])
    Alfa=float(paramValues['ALFA'])
    zlow=float(paramValues['ZLOW'])
    zup=float(paramValues['ZUP'])
    zerr=float(paramValues['ZERR'])
    zderr=float(paramValues['ZDERR'])
    mlim=float(paramValues['MLIM'])
    merr=float(paramValues['MERR'])
    mderr=float(paramValues['MDERR'])
    color=float(paramValues['C'])
    colord=float(paramValues['CD'])
    cerr=float(paramValues['CERR'])
    cderr=float(paramValues['CDERR'])
    density=float(paramValues['DENSITY'])
    nobjects=[float(i) for i in paramValues['NOBJECTS']]

    for numObjects in nobjects:
        if verbose: output.write('Dentro del bucle. '+str(datetime.datetime.today())+'\n')

        area=numObjects/density
       
        #Names of files
        simulationName="simlf_"+tipo+"_mstar"+str(Mstar)+"_phistar"
        simulationName+=str(Phistar)+"_alfa"+str(Alfa)+"_zlow"+str(zlow)+"_zup"+str(zup)
        simulationName+="_zerr"+str(zerr)+"_zderr"+str(zderr)
        simulationName+="_mlim"+str(mlim)+"_merr"+str(merr)+"_mderr"+str(mderr)
        simulationName+="_c"+str(color)+"_cd"
        simulationName+=str(colord)+"_cerr"+str(cerr)
        simulationName+="_cderr"+str(cderr)
        simulationName+="_area"+str(area)+"/"
   
        catDir=simulationsDir+"/SimulatedCatalogs/"+simulationName
        simulatedCatalogFile=catDir+str(simulationNumber)+".cat"

        lfDir=simulationsDir+'/ComputedLF/STY_wC/'+simulationName
        if not(os.access(lfDir,os.F_OK)):
            os.makedirs(lfDir)
            if verbose: output.write("Creating "+lfDir+"\n")
        computedSTY_wCFile=lfDir+str(simulationNumber)+".lf"

        # done y lock files
        doneDir=simulationsDir+"/LockFiles/STY_wC/"+simulationName
        if not(os.access(doneDir,os.F_OK)):
            os.makedirs(doneDir)
            if verbose: output.write("Creating "+doneDir+"\n")

        doneCatFileName =simulationsDir+"/LockFiles/Cat/"+simulationName
        doneCatFileName+=str(simulationNumber)+".done"
        doneSTY_wCFileName=doneDir+str(simulationNumber)+".done"
        lockSTY_wCFileName=doneDir+str(simulationNumber)+".lock"

        comando ="LF_computeSTY_wC_Mag.csh "+simulatedCatalogFile+" "
        comando+=str(color)+" "+str(colord)+" "+str(mlim)
        comando+=" "+str(zlow)+" "+str(zup)+" "+str(area)
	comando+=" "+computedSTY_wCFile

        doneFile=0
        doneFile=os.access(doneCatFileName,os.F_OK)

        lockFileName=lockSTY_wCFileName # para el capturador de sigterm
        global PID
        PID=False

        if doneFile:
            # mandamos el catálogo
            estatusSTY_wC=secureCompute('STY_wC',doneSTY_wCFileName,lockSTY_wCFileName,comando,computedSTY_wCFile)
            if estatusSTY_wC==1:
                if verbose: output.write("Problemas con STY_wC\n")
	        continue # seguimos con otro
            else:
                if verbose: output.write("Todo correcto con STY_wC\n")
	 
            if verbose: output.write("Ya está hecho STY_wC para nobjects ="+str(numObjects)+"\n")
            continue
        else:
            if verbose: output.write("No se puede hacer STY_wC porque no está el catálogo.\n")

    if verbose: output.write("Iter ="+str(simulationNumber)+"\n")

    # Cuando termina el bucle en número de galaxias
    if verbose: output.write("\nLF_launch_single_STY_wC.py termina aquí\n")
    if verbose: output.write("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
    output.close()
    sys.exit(0)
   
if __name__ == "__main__":
    main()
