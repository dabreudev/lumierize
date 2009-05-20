#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# TITLE: Script para analizar las simulaciones con LumFunc_GOYA
# CREATED: 20080325 dabreu
# Este script analiza el resultado de una simulaci�n.
#
# MODIFIED: 20080402 dabreu
# Added plots of ML_MAX vs M_STAR
#
# MODIFIED: 20080425 dabreu
# Added F (Fisher)
#
# Uso: LF_analyze_simulations.py -m VVmax simulationParameterFile
#

"""Este script analiza el resultado de una simulaci�n"""

import sys
from optparse import OptionParser
import os
import catg
from scipy import stats
import numpy
import pylab
from utilidades import sigmaDelCuartil
from math import sqrt
mean=stats.mean
median=stats.median
g=pylab

simulationsDir= os.getenv("LF_SIM")

#Funci�n para leer los par�metros de la simulaci�n
def readParams(file):
    """readParams(file) -> devuelve un diccionario con los datos del fichero de
                        par�metros en un diccionario y los valores como
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
	    if len(trozos[1:])>1: #par�metro con m�s de un valor
	        parametros[trozos[0]]=trozos[1:] #se devuelve como una lista
	    else:
	        parametros[trozos[0]]=trozos[1]

    return parametros

def main():
    usage = "usage: %prog [options] simulationParameterFile \nDefault values: out=stdout, verbose=1, method=VVmax, plots=no"

    parser = OptionParser(usage)
    parser.add_option("-o", "--out", dest="outfilename", default="stdout",
                      help="Output to FILE [default: %default]", metavar="FILE")
    parser.add_option("-v", "--verbose", dest="verbose", default=1,
                      help="Set verbose level to INT [default: %default", metavar="INT")
    parser.add_option("-m", "--method", dest="type", default='VVmax',
                choices=["VVmax","STY","STY_wC","STY_wC_errColor","STY_errMag"],
                help="Simulation TYPE [default: %default]", metavar="TYPE")
    parser.add_option("-p", "--plots", dest="plots", default='no',
                      choices=["no","ps","png"],
                      help="Plots of TYPE [default: %default]", metavar="TYPE")
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
    if options.plots=="no":
        plotType=False
    else:
        plotType=options.plots
   
    paramValues=readParams(simulationParameterFile)

    niters=int(paramValues['NITER'])   
    tipo=paramValues['TYPE']
    H0=float(paramValues['H0'])
    OM=float(paramValues['OM'])
    OL=float(paramValues['OL'])
    Mstar0=float(paramValues['MSTAR'])
    Phistar0=float(paramValues['PHISTAR'])
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
    density=float(paramValues['DENSITY'])
    nobjects=[float(i) for i in paramValues['NOBJECTS']]

    valores=['length','nonan','mean','median','sigcuartil']
    parametros=['Mstar','Phistar','Alpha']
    parametros0={'Alpha0':Alfa0,'Phistar0':Phistar0,'Mstar0':Mstar0}
    llaves={'Mstar':'M_STAR','Phistar':'PHISTAR','Alpha':'ALPHA',
            'ErrMstar':'ERR_M_STAR','ErrPhistar':'ERR_PHISTAR',
  	    'ErrAlpha':'ERR_ALPHA'}

    # variables diccionario con los resultados
    resultados={}
    MLmaxData={}
    MLmaxData['MLmax']=[]
    histogramas={}
    for parametro in parametros:
        for valor in valores:
            for i in ['','Err']:
	        resultados[valor+i+parametro]=[]
    for parametro in parametros:
        MLmaxData[parametro]=[]
        for i in ['','Err']:
            histogramas[i+parametro]=[]         

    for numObjects in nobjects:
        if verbose: output.write('Dentro del bucle. '+str(int(numObjects))+'\n')
        area=numObjects/density

        #Names of files
        simulationName="_H0"+str(H0)+"_OM"+str(OM)+"_OL"+str(OL)
        simulationName+="_mstar"+str(Mstar0)+"_phistar"+str(Phistar0)
        simulationName+="_alfa"+str(Alfa0)
        simulationName+="_zlow"+str(zlow)+"_zup"+str(zup)
        simulationName+="_zerr"+str(zerr)+"_zderr"+str(zderr)
        simulationName+="_mlim"+str(mlim)+"_deltm"+str(deltamag)       
        simulationName+="_merr"+str(merr)+"_mderr"+str(mderr)
        simulationName+="_c"+str(color)+"_cd"+str(colord)
        simulationName+="_cerr"+str(cerr)+"_cderr"+str(cderr)
        simulationName+="_area"+str(area)+"/"


        resultsDir=simulationsDir+"/Results/"+method+"/cosecha"+simulationName+"/"

        dataFile=resultsDir+method+"_numObjects_"+str(int(numObjects))+".cat"

        if verbose: output.write("Reading data from "+dataFile+"\n")
        data=catg.rcat(dataFile)
        if verbose: output.write("done.\n")

        MLmaxData['MLmax'].append(data['ML_MAX'])
        for parametro in parametros:
            MLmaxData[parametro].append(data[llaves[parametro]])
            for i in ['','Err']:
	        j=i+parametro
                resultados['length'+j].append(len(data[llaves[j]]))
                # para evitar los nan
                d=[float(k) for k in data[llaves[j]]]
                d=numpy.asarray(d)
                l=numpy.where(d <= numpy.max(d))
                resultados['nonan'+j].append(len(l[0]))
                try:
                    resultados['mean'+j].append(mean(data[llaves[j]]))
                except:
                    resultados['mean'+j].append(float('nan'))
                try:
                    resultados['median'+j].append(median(data[llaves[j]]))
                except:
                    resultados['median'+j].append(float('nan'))
                try:
                    resultados['sigcuartil'+j].append(sigmaDelCuartil(data[llaves[j]]))
                except:
                    resultados['sigcuartil'+j].append(float('nan'))

                try:
                    histo=pylab.histogram(data[llaves[j]],bins=30)#,defaultlimits=limits)
                except:
                    histo=[[0],[0]] 
                histogramas[j].append(histo)

    #sys.exit()

    #resCat=catg.newcat('N_OBJECTS', nobjects) # to put results in a catalogue

    totalT=[]

    if verbose: output.write("Output served in twiki format for your comfort.\n\n")   
    output.write("---++Cosecha"+simulationName+"\n\n")
    for i in range(len(nobjects)):

        # Pintamos resultados
        output.write("---+++Results for numObjects = "+str(int(nobjects[i]))+"\n\n")
        output.write("   * Number of simulations: "+str(resultados["lengthMstar"][i])+"\n")
        output.write("   * nonan -> number of values != NaN\n")
        output.write("   * T -> (median-orig)/(sigmaDelCuartil/sqrt(N))\n\n")

        output.write("| Parameter | nonan | nonan(Err) | mean | median | sigma | median(Err) | median-orig | T |\n")

        format="%2.4f"
        
        partialT=[]

        for parametro in parametros:
            diff=resultados["median"+parametro][i]-parametros0[parametro+'0']
            sigma=resultados["sigcuartil"+parametro][i]
            output.write("| "+parametro)
            output.write(" | %i" % resultados["nonan"+parametro][i])
            output.write(" | %i" % resultados["nonanErr"+parametro][i])
            output.write(" | "+format % resultados["mean"+parametro][i])
            output.write(" | "+format % resultados["median"+parametro][i])
            output.write(" | "+format % sigma)
            #print "kk "+str(resultados["medianErr"+parametro][i])+"\n"
            try:
                output.write(" | "+format % float(resultados["medianErr"+parametro][i]))
            except:
                output.write(" | "+str(resultados["medianErr"+parametro][i]))
            output.write(" | "+format % diff)
            t=diff*sqrt(resultados["length"+parametro][i])/sigma
            output.write(" | "+format % t)
            output.write(" |\n")
            partialT.append(abs(t))
        output.write("\n")
        totalT.append(mean(partialT))

    output.write("---+++ Mean value for T.\n\n")
    output.write("| numobjects | mean(abs(T)) |\n")
    for i in range(len(nobjects)):
        output.write("| "+str(int(nobjects[i])))
        output.write(" | "+format % totalT[i])
        output.write(" |\n")
    output.write("\n")

    if plotType:
        if verbose: output.write("You have selected plots in "+plotType+" format.\n")

        # ML_MAX vs M_STAR, ...
        if verbose: output.write("\nMLmax plots.\n============\n")
        for i in range(len(nobjects)):
            ni=str(int(nobjects[i]))
            if verbose: output.write("\nnumObjects = "+ni+"\n")
            if verbose: output.write("------------\n")
            for parametro in parametros:
                #print MLmaxData['MLmax'][i]
                #print MLmaxData[parametro][i]
                graphName =resultsDir+"MLmax_"+parametro+"_numObjects_"+ni
                graphName+="."+plotType
                t ="For numObjects "+ni+" ("
                t+=str(resultados["length"+parametro][i])+")"
                pylab.plot(MLmaxData[parametro][i],MLmaxData['MLmax'][i],'.')
                pylab.title(t)
                pylab.xlabel(parametro)
                pylab.ylabel('MLmax')
                pylab.savefig(graphName)
                pylab.clf()
                output.write("Graph of MLmax "+parametro+" in: "+graphName+"\n")

        # Histograms
        if verbose: output.write("\nHistograms.\n===========\n")
        for i in range(len(nobjects)):
            ni=str(int(nobjects[i]))
            if verbose: output.write("\nnumObjects = "+ni+"\n")
            if verbose: output.write("------------\n")
            for parametro in parametros:
                for j in ['','Err']:
                    graphName=resultsDir+"Histo_"+j+parametro+"_numObjects_"+ni+"."+plotType
                    t ="Histogram of "+parametro+j+" for numObjects "+ni+" ("
                    t+=str(resultados["length"+j+parametro][i])+" simulations)"
                    h=histogramas[j+parametro][i]
                    try:
                        pylab.bar(h[1],h[0],width=0.9*(h[1][1]-h[1][0]),ec='blue',fill=True)
                        pylab.title(t)
                        pylab.xlabel(j+parametro)
                        pylab.savefig(graphName)
                        pylab.clf()
                        output.write("Graph of "+j+parametro+" in: "+graphName+"\n")
                    except:
                        output.write("I could not do graph of "+j+parametro+"\n")

        if verbose: output.write("A ration of histograms ready for you!\n")

    if verbose: output.write("All done. I hope you are happier now.\n")
    output.close()
    sys.exit()
   
    #output.write('MStar0: '+str(MStar0)+'\n')
    #for valor in valores:
    #    for i in ['','Err']:
    #	    j=valor+i+'MStar'
    #        output.write(j+': '+str(resultados[j])+'\n')
    #        catg.addcol(resCat, j, resultados[j])
    #
    #output.write('PhiStar0: '+str(PhiStar0)+'\n')
    #for valor in valores:
    #    for i in ['','Err']:
    #	    j=valor+i+'PhiStar'
    #        output.write(j+': '+str(resultados[j])+'\n')
    #        catg.addcol(resCat, j, resultados[j])

    #output.write('Alfa0: '+str(Alfa0)+'\n')
    #for valor in valores:
    #    for i in ['','Err']:
    #	    j=valor+i+'Alpha'
    #        output.write(j+': '+str(resultados[j])+'\n')
    #        catg.addcol(resCat, j, resultados[j])
    #
    #catg.writecat(resCat, 'resultados.cat', overwrite=1)
    #output.write('Resultados en: "resultados.cat"\n')
      
    # histogramas en cada par�metro entre +- 4 � 5 sigmas del
    # cuartil (la resta)
   
    #for parametro in parametros:
    #    for i in ['','Err']:
    #	    j=0
    #        for data in histogramas[i+parametro+metodo]:
    #	        numObjects=nobjects[j]
    #	        yHisto=data[0]
    #	        xHisto=numpy.asarray(range(len(yHisto)))*data[2]+data[1]
    #           plot=g.bar(xHisto, yHisto, width=data[2])
    #	        grafFileName='Results/'+i+parametro+metodo+str(numObjects)+'.ps'
    #	        g.title(i+parametro+' '+metodo+' '+str(numObjects)+' galaxias/simulacion')
    #	        output.write('Gr�fica en: '+grafFileName+'\n')
    #           g.savefig(grafFileName)
    #	        g.savefig(grafFileName[:-5]+'.png')
    #            g.clf()
    #	        j+=1

    sys.exit()
      
if __name__ == "__main__":
    main()
