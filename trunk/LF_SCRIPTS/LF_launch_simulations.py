#! /usr/bin/env python
# -*- coding: ISO-8859-1 -*-

# TITLE: Script en python para realizar simulaciones
# CREATED: 20060913 dabreu
# Pequeño script para realizar simulaciones que se encarga de llamar N veces un
# comando y manejar los SIGTERM (muy útil en el caso de condor)
#
# MODIFIED: 20061004 dabreu
# Se añade el TIMEOUT al pexpect para capturar un proceso que tarde demasiado y
# matarlo correctamente
#
# MODIFIED: 20061212 dabreu
# Cambios en los parámetros
#
# MODIFIED: 20080323 dabreu
# Modificado para que sea genérico y acepte otro argumento.
#
# Uso: LF_launch_simulations.py -p parFile nIters
#
"""Script para realizar simulaciones"""

import sys
from optparse import OptionParser
import os
import pexpect
import signal
import datetime
import time

def main():

   def signal_handler(signal, frame):
      """signal_handler(signal, frame) -> Para poder capturar el sigterm matar el
                                       el proceso con SIGTERM

         Para que se produzca la captura, la función debe ser llamada con:
         signal.signal(signal.SIGTERM, signal_handler)"""
      output.write("Recibí un SIGTERM (habló launch_simulations)\n")
      p.kill(15)
      time.sleep(0.1) # to wait for the child process
      sys.exit(-15)

   signal.signal(signal.SIGTERM, signal_handler)

   # Aquí empieza main
   usage = "usage: %prog [options] singleSimulatorName nIter \nDefault values: out=stdout, verbose=1 par=default.par"
   parser = OptionParser(usage)
   parser.add_option("-o", "--out", dest="outfilename",
                     help="Output to FILE", metavar="FILE")
   parser.add_option("-v", "--verbose", dest="verbose",
                     help="Set verbose level to INT", metavar="INT")
   parser.add_option("-p", "--par", dest="param",
                     help="Parameter filename (from WORK dir)", metavar="FILE")
   parser.set_default("outfilename","stdout")
   parser.set_default("param","default.par")
   parser.set_default("verbose",1)
   
   (options, args) = parser.parse_args()
   
   if len(args)<2:
      print "Error in args. Type --help for help\n"
      sys.exit(1)

   scriptName=args[0]
   nIter=int(args[1])
   parFile=options.param
   verbose=options.verbose

   if options.outfilename=="stdout":
      output=sys.stdout
   else:
      output=open(options.outfilename, "w")
   
   for iter in range(nIter):
      if verbose: output.write('Dentro del bucle. '+str(datetime.datetime.today())+'\n')
      comando=scriptName+" -p "+parFile+" "+str(iter)
      if verbose: output.write("Antes del pexpect.spawn\n")
      p=pexpect.spawn(comando, timeout=21600) # 6 horas como máximo!
      
      queHaPasado=p.expect([pexpect.EOF,pexpect.TIMEOUT])
      
      if queHaPasado==0:
         estatus=os.getenv("status")
         if verbose: output.write(p.before)
         if verbose: output.write("Después del p.expect\n")
      if queHaPasado==1:
         if verbose: output.write("He obtenido un TIMEOUT\n")
	 if verbose: output.write(p.before)
         p.kill(15)

      # hay que matarlo y hacer unos sleep, porque si no aparecen
      # como <defunct>
      output.write("Cerramos proceso "+str(p.pid)+".\n")
      p.close()
      time.sleep(0.1)
      p.kill(15)
      time.sleep(0.1)
      p.kill(9)

   output.write("\nTodo correcto\n")    
   output.close()
   sys.exit()

if __name__ == "__main__":
   main()
