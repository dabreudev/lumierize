#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# TITLE: Script para borrar ficheros de una simulación
# CREATED: 20060918 dabreu
# Este script borra los ficheros que se obtienen como resultado de una
# simulación. Es necesario debido a que el número de ficheros es muy grande y se
# obtiene un error en "rm".
#
# MODIFIED: 20080324 dabreu
# Added new parameters.
#
# Uso: LF_clean.py -p paramFileName todos
#

"""Script para borrar ficheros de una simulación"""

import sys
from optparse import OptionParser
import os

parentDir=os.getenv("LF_SIM")

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
   usage = "usage: %prog [options] string \nDefault values: par=None exclude=None out=stdout verbose=1"
   usage+="\n\n"
   usage+="Poniendo 'todos' como agumento se asume que se borran todos los archivos\n"
   usage+="pertenecientes a la simulación definida por el fichero de parámetros"
   parser = OptionParser(usage)
   parser.add_option("-o", "--out", dest="outfilename",
                     help="Output to FILE", metavar="FILE")
   parser.add_option("-v", "--verbose", dest="verbose",
                     help="Set verbose level to INT", metavar="INT")
   parser.add_option("-p", "--par", dest="param",
                     help="Parameter filename (from WORK dir)", metavar="FILE")
   parser.add_option("-e", "--exclude", dest="exclude",
                     help="Excluded directory (from WORK dir)", metavar="DIR")
   parser.set_default("outfilename","stdout")		     
   parser.set_default("param",None)
   parser.set_default("exclude",None)
   parser.set_default("verbose",1)
   
   (options, args) = parser.parse_args()

   verbose=options.verbose
   parameterFileName=options.param
   excludedDir=options.exclude

   if len(args)<1:
      print "Type --help for help\n"
      sys.exit()
   else:
      string=args[0]
   
   if options.outfilename=="stdout":
      output=sys.stdout
   else:
      output=open(options.outfilename, "w")

   #directorios=['ComputedLF','LockFiles','Logs','SimulatedCatalogs']
   directorios=['ComputedLF','LockFiles','SimulatedCatalogs']
   
   if excludedDir != None:
      try:
         directorios.remove(excludedDir)
      except:
         pass
   
   if parameterFileName == None:
      if string == 'todos':
         # Se borran TODOS los ficheros
         for directorio in directorios:
            nDel=0
            files=os.listdir(parentDir+'/'+directorio+'/')
            for file in files:
               nDel+=1
               os.remove(parentDir+'/'+directorio+'/'+file)

            if(verbose): output.write("Borrados "+str(nDel)+" los ficheros de "+directorio+"\n")
         
      else:
         # Se borran todos los ficheros que contengan string
	 for directorio in directorios:
            nDel=0
	    files=os.listdir(parentDir+'/'+directorio+'/')
	    for file in files:
	       if string in file:
                  os.remove(parentDir+'/'+directorio+'/'+file)
                  nDel+=1

            if(verbose): output.write("Borrados "+str(nDel)+" ficheros de "+directorio+"\n")

   else:
      # habrá que leer entonces el fichero de parámetros
      paramValues=readParams(parameterFileName)

      #tipo=paramValues['TYPE']
      H0=float(paramValues['H0'])
      OM=float(paramValues['OM'])
      OL=float(paramValues['OL'])
      MStar0=float(paramValues['MSTAR'])
      PhiStar0=float(paramValues['PHISTAR'])
      Alfa0=float(paramValues['ALFA'])
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
      #density=float(paramValues['DENSITY'])
      #nobjects=[float(i) for i in paramValues['NOBJECTS']]

      #Names of files
      simulationName+="_H0"+str(H0)+"_OM"+str(OM)+"_OL"+str(OL)
      simulationName+="_mstar"+str(Mstar)+"_phistar"+str(Phistar)                              simulationName+="_alfa"+str(Alfa)
      simulationName+="_zlow"+str(zlow)+"_zup"+str(zup)
      simulationName+="_zerr"+str(zerr)+"_zderr"+str(zderr)
      simulationName+="_mlim"+str(mlim)+"_deltm"+str(deltamag)                                 simulationName+="_merr"+str(merr)+"_mderr"+str(mderr)
      simulationName+="_c"+str(color)+"_cd"+str(colord)
      simulationName+="_cerr"+str(cerr)+"_cderr"+str(cderr)
      simulationName+="_area"+str(area)+"/"

      if string == 'todos':
         # Se borran TODOS los ficheros de esa simulación
         for directorio in directorios:
            nDel=0
            files=os.listdir(parentDir+'/'+directorio+'/')
            for file in files:
               if simulationName in file:
                  nDel+=1
                  os.remove(parentDir+'/'+directorio+'/'+file)
                  if verbose > 2: ouput.write("Borrado: "+directorio+'/'+file+'\n')

            if(verbose): output.write("Borrados "+str(nDel)+" ficheros de "+directorio+"\n")
               
      else:
         # Se borran todos los ficheros de esa simulación que contengan string
	 for directorio in directorios:
            nDel=0
	    files=os.listdir(parentDir+'/'+directorio+'/')
	    for file in files:
	       if simulationName in file and string in file:
                  nDel+=1
	          os.remove(parentDir+'/'+directorio+'/'+file)
                  if verbose > 2: ouput.write("Borrado: "+directorio+'/'+file+'\n')

            if(verbose): output.write("Borrados "+str(nDel)+" ficheros de "+directorio+"\n")

   if(verbose): output.write("Se ha terminado de hacer el clean de '"+string+"'\n")
   sys.exit()

if __name__ == "__main__":
   main()
