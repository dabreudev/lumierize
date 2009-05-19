#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# TITLE: Script to check the status of a simulation
# CREATED: 20080324 dabreu
# This script reads all output files and gives the status of a simulation of a
# given parameter file
#
# MODIFIED: 20080422 dabreu
# LockFiles structure modified
#  
# Use: LF_status.py simulationParameterFile
#

"""Script to check the status of a simulation"""

import sys
from optparse import OptionParser
import os
import signal
import time
import commands
import datetime

estatus=0
simulationsDir=os.getenv("LF_SIM")

#To read parameter of a simulation
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

def giveMeNumber(data,name):
    """giveMeNumber(data,name) -> count number of name in data[0]"""

    number=0

    dir=data[0]+name
    files=os.listdir(dir)
    for file in files:
        if data[1] in file:
            number+=1

    return number

def main():

    # Aquí empieza main
    usage = "usage: %prog [options] simulationParameterFile\nDefault values: out=stdout, verbose=1"
    parser = OptionParser(usage)
    parser.add_option("-o", "--out", dest="outfilename",
                     help="Output to FILE", metavar="FILE")
    parser.add_option("-v", "--verbose", dest="verbose",
                     help="Set verbose level to INT", metavar="INT")
    parser.add_option("-l", "--locks", dest="process", default='0',
                     choices=['0','1'],
                     help="Look only for process  [default: %default]",
                     metavar="INT")
    parser.add_option("-m", "--method", dest="method", default="all",
                     choices=['all','Cat','VVmax','STY','STY_wC','STY_errMag'],
                     help="Show only method METHOD [default: %default]",
                     metavar="METHOD")
    parser.set_default("outfilename","stdout")
    parser.set_default("param","default.par")
    parser.set_default("verbose",1)
   
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Error in args. Type --help for help\n"
        sys.exit(1)

    simulationParameterFile=args[0]
    verbose=options.verbose
    methodSelected=options.method   
    if options.process=='0':
        process=False
    else:
        process=True

    if options.outfilename=="stdout":
        output=sys.stdout
    else:
        output=open(options.outfilename, "w")
      
    parameterFile=simulationsDir+"/"+simulationParameterFile
    paramValues=readParams(parameterFile)

    niter=int(paramValues['NITER'])   
    tipo=paramValues['TYPE'] # "m" or "l", but we are using only "m"
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

    nameList=[]

    for numObjects in nobjects:
        area=numObjects/density
       
        #Names of files
        simulationName="simlf_"+tipo
        simulationName+="_H0"+str(H0)+"_OM"+str(OM)+"_OL"+str(OL)
        simulationName+="_mstar"+str(Mstar)+"_phistar"+str(Phistar)
        simulationName+="_alfa"+str(Alfa)
        simulationName+="_zlow"+str(zlow)+"_zup"+str(zup)
        simulationName+="_zerr"+str(zerr)+"_zderr"+str(zderr)
        simulationName+="_mlim"+str(mlim)+"_deltm"+str(deltamag)
        simulationName+="_merr"+str(merr)+"_mderr"+str(mderr)
        simulationName+="_c"+str(color)+"_cd"+str(colord)
        simulationName+="_cerr"+str(cerr)+"_cderr"+str(cderr)
        simulationName+="_area"+str(area)+"/"

        nameList.append(simulationName)

    if methodSelected=='all':
        methods=['Cat','VVmax','STY','STY_wC','STY_errMag','STY_wC_errColor']
    else:
        methods=[methodSelected]

    if process:
        allLocks=simulationsDir+'/LockFiles/'
        methodsData={'Cat':[allLocks+'Cat/','.lock'],
                     'VVmax':[allLocks+'VVmax/','.lock'],
                     'STY':[allLocks+'STY/','.lock'],
                     'STY_wC':[allLocks+'STY_wC/','.lock'],
                     'STY_wC_errColor':[allLocks+'STY_wC_errColor/','.lock'],
                     'STY_errMag':[allLocks+'STY_errMag/','.lock']}
    else:
       	allCatalogs=simulationsDir+'/SimulatedCatalogs/'
       	allLFs=simulationsDir+'/ComputedLF/'
       	methodsData={'Cat':[allCatalogs,''],'VVmax':[allLFs+'VVmax/','.lf'],
              	      'STY':[allLFs+'STY/','.lf'],
              	      'STY_wC':[allLFs+'STY_wC/','.lf'],
              	      'STY_wC_errColor':[allLFs+'STY_wC_errColor/','.lf'],
              	      'STY_errMag':[allLFs+'STY_errMag/','.lf']}

    output.write("Status of simulation "+simulationParameterFile+"\n")
    output.write("==============================================\n\n")

    for method in methods:
        output.write(method+"\n")
        output.write("------------------\n")
        totalNumber=0
        for i in range(len(nobjects)):
            output.write("Ngals = "+str(nobjects[i])+"\t-> ")
            number=giveMeNumber(methodsData[method],nameList[i])
            totalNumber+=number
            if process:
                output.write(str(number)+"\trunning\n")
            else:
                output.write(str(number)+"\t/ "+str(niter)+"\tdone\n")
        output.write("\nTotal\t\t->"+str(totalNumber)+"\n\n")

    output.close()
    sys.exit(0)
   
if __name__ == "__main__":
    main()
