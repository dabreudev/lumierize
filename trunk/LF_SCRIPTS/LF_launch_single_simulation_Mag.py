#! /usr/bin/env python

# TITLE: Script para lanzar una sola simulacion con LumFunc_GOYA de STY_Mag
# CREATED: 20060907 dabreu
# Este script lanza una simulación con unos parámetros dados, calcula
# la LF por STY y por VVmax.
#
# MODIFIED: 20061212 dabreu
# Ajuste de cambios en los parámetros
#
# Uso: LF_launch_single_simulation_Mag.py -p simulationParameterFile numberOfIter
#

"""Script para lanzar una sola simulacion con LumFunc_GOYA"""

import sys
from optparse import OptionParser
import os
import signal
import time
import commands
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

      output.write("Recibí un SIGTERM (habló launch_single_simulation)\n")
      try:
         os.remove(lockFileName)
         output.write("El lockfile se borró correctamente\n")
      except:
         output.write("La petó al intentar borrar el lockfile\n")
         pass
      
      sys.exit(1)

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
   	    (estatus, salida)=commands.getstatusoutput(comando)
   	    if verbose: output.write(salida)
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
   usage = "usage: %prog [options] nIter \nDefault values: out=stdout, par=default.par, verbose=1"
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
      
   parameterFile=simulationsDir+simulationParameterFile
   paramValues=readParams(parameterFile)

   niter=int(paramValues['NITER'])   
   tipo=paramValues['TYPE']
   Mstar=float(paramValues['MSTAR'])
   Phistar=float(paramValues['PHISTAR'])
   Alfa=float(paramValues['ALFA'])
   zlow=float(paramValues['ZLOW'])
   zup=float(paramValues['ZUP'])
   zerr=float(paramValues['ZERR'])
   zderr=float(paramValues['ZDERR'])
   mlim=float(paramValues['MLIM'])
   merr=float(paramValues['MERR'])
   errdev=float(paramValues['ERRDEV'])
   density=float(paramValues['DENSITY'])
   nobjects=[float(i) for i in paramValues['NOBJECTS']]

   for numObjects in nobjects:
      if verbose: output.write('Dentro del bucle. '+str(datetime.datetime.today())+'\n')

      area=numObjects/density
       
      #Names of files
            
      simulationName="_"+tipo+"_"+str(simulationNumber)+"_mstar"+str(Mstar)+"_phistar"
      simulationName+=str(Phistar)+"_alfa"+str(Alfa)+"_zlow"+str(zlow)+"_zup"+str(zup)
      simulationName+="_zerr"+str(zerr)+"_zderr"+str(zderr)
      simulationName+="_mlim"+str(mlim)+"_merr"+str(merr)+"_errdev"+str(errdev)
      simulationName+="_area"+str(area)
   
      simulatedCatalogFile=simulationsDir+"SimulatedCatalogs/simlf"+simulationName+".cat"
      computedSTY_p_MFile=simulationsDir+"ComputedLF/simlf_STY_p_M"+simulationName+".lf"
      computedVVmaxFile=simulationsDir+"ComputedLF/simlf_VVmax"+simulationName+".lf"

      # done y lock files
      doneFileName=simulationsDir+"LockFiles/simlf"+simulationName+"All.done"
      doneFile=0
      
      doneFile=os.access(doneFileName,os.F_OK)
      if not(doneFile): # El paso no está hecho

         doneCatFileName=simulationsDir+"LockFiles/simlf"+simulationName+"Cat.done"
         lockCatFileName=simulationsDir+"LockFiles/simlf"+simulationName+"Cat.lock"
         doneVVmaxFileName=simulationsDir+"LockFiles/simlf"+simulationName+"VVmax.done"
         lockVVmaxFileName=simulationsDir+"LockFiles/simlf"+simulationName+"VVmax.lock"
         doneSTYFileName=simulationsDir+"LockFiles/simlf"+simulationName+"STY.done"
         lockSTYFileName=simulationsDir+"LockFiles/simlf"+simulationName+"STY.lock"

         comandoCat ="LF_generate_catalog_Mag.csh "+str(Mstar)+" "+str(Phistar)
         comandoCat+=" "+str(Alfa)+" "+str(zlow)+" "+str(zup)+" "+str(zerr)
	 comandoCat+=" "+str(zderr)+" "+str(mlim)+" "
         comandoCat+=str(merr)+" "+str(errdev)+" "+str(area)+" "
         comandoCat+=simulatedCatalogFile

         comandoVVmax ="LF_computeVVmax_Mag.csh "+simulatedCatalogFile+" "
         comandoVVmax+=str(zlow)+" "+str(zup)+" "+str(area)+" "+str(mlim)
         comandoVVmax+=" "+computedVVmaxFile

         comandoSTY ="LF_computeSTY_Mag.csh "+simulatedCatalogFile+" "
         comandoSTY+=str(mlim)+" "+str(zlow)+" "+str(zup)+" "+str(area)
         comandoSTY+=" "+computedSTY_p_MFile
		  
         # mandámos el catálogo
         estatusCat=secureCompute('catálogo',doneCatFileName,lockCatFileName,comandoCat,simulatedCatalogFile)
         if estatusCat==1:
            if verbose: output.write("Problemas con el catálogo\n")
	    continue # si no hay catálogo no se puede hacer vvmax ni sty
         else:
            if verbose: output.write("Todo correcto con el catálogo\n")
	 
         # mandamos VVMAX
         estatusVVmax=secureCompute('VVmax',doneVVmaxFileName,lockVVmaxFileName,comandoVVmax,computedVVmaxFile)
         if estatusVVmax==1:
            if verbose: output.write("Problemas con vvmax\n")
         else:
            if verbose: output.write("Todo correcto con vvmax\n")

         # mandamos STY
         estatusSTY=secureCompute('STY',doneSTYFileName,lockSTYFileName,comandoSTY,computedSTY_p_MFile)
         if estatusSTY==1:
            if verbose: output.write("Problemas con STY\n")
         else:
            if verbose: output.write("Todo correcto con STY\n")

         # Hacemos el done si todo terminó correcto
         if estatusVVmax==0 and estatusSTY ==0:
	    if verbose: output.write("Ponemos el done\n")
	    os.system('touch '+doneFileName)

      else: # Ya está todo hecho
         if verbose: output.write("Ya está hecho para nobjects ="+str(numObjects)+"\n")
         continue

   # Cuando termina el bucle en número de galaxias
   if verbose: output.write("\nLF_launch_single_simulations_Mag.py termina aquí\n")
   if verbose: output.write("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
   output.close()
   sys.exit(0)
   
if __name__ == "__main__":
   main()
