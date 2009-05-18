#! /bin/csh

# TITLE: Script para generar catálogos con LumFunc_GOYA wC
# CREATED: 20080212 dabreu
# Pequeño script para generar un catálogo en magnitudes para realizar
# simulaciones utilizando la opción I de LumFunc_GOYA
#
# MODIFIED: 20090518 dabreu
# Support for new cosmology
#
# Argumentos [valor por defecto]:
#  - H0 [70]
#  - Omega matter [0.3]
#  - Omega lambda [0.7]
#  - Mstar [-20.4]
#  - Phistar [0.0033]
#  - Alfa [-1.3]
#  - Z up [0]
#  - Z low [0]
#  - z error medio [0]
#  - z desviación del error [0]
#  - magnitud limite [25] (salen ~2800 objetos)
#  - deltamag [0.0]
#  - error medio [0.0]
#  - desviación del error [0.0]
#  - mean Color [0.0]
#  - stddev Color [0.0]
#  - mean errorColor [0.0]
#  - stddev errorColor [0.0]
#  - área del survey [0.1]
#  - fichero de salida [kk.cat]
#
# Solo hay dos formas de llamarlo correctamente, dando todos los argumentos o
# ninguno, en cuyo caso tomarán los valores por defecto
#

set estatus=0 #para guardar el status de salida de LumFunc_GOYA

if ($1 == "--help") then
   echo "Uso: LF_generate_catalog_Mag_wC.csh H0 OM OL Mstar Phistar Alfa zlow zup zError zdError magLim deltaMag mError dError meanColor stddevColor meanErrorColor meanStddevColor area outFile"
   echo ""
   echo "Example: LF_generate_catalog_Mag_wC.csh 70 0.3 0.7 -20.4 0.0033 -1.3 0 0 0 0 25 0 0 0 0 0 0 0 0.1 kk.cat"
   exit
endif

if ($1 == "") then
   set H0=70
   set OM=0.3
   set OL=0.7
   set Mstar=-20.4
   set Phistar=0.0033
   set Alfa=-1.3
   set zlow=0
   set zup=0
   set zError=0
   set zdError=0
   set magLim=25
   set deltaMag=0.0
   set mError=0.0
   set dError=0.0
   set meanColor=0.0
   set stddevColor=0.0
   set meanErrorColor=0.0
   set stddevErrorColor=0.0
   set area=0.1
   set outFile="kk.cat"
else
   set H0=$1
   set OM=$2
   set OL=$3
   set Mstar=$4
   set Phistar=$5
   set Alfa=$6
   set zlow=$7
   set zup=$8
   set zError=$9
   set zdError=$10
   set magLim=$11
   set deltaMag=$12
   set mError=$13
   set dError=$14
   set meanColor=$15
   set stddevColor=$16
   set meanErrorColor=$17
   set stddevErrorColor=$18
   set area=$19
   set outFile=$20
endif

$OPERA_INST/LumFunc_GOYA << COMANDOS
I
$H0
$OM
$OL
$Mstar
$Phistar
$Alfa
$zlow
$zup
$zError
$zdError
$magLim
$deltaMag
10
$mError
$dError
$meanColor
$stddevColor
$meanErrorColor
$stddevErrorColor
$area
0
$outFile
e
COMANDOS

set estatus=$status

echo ""
echo "-------------------------"
echo $0 $*":"
echo "H0: $H0"
echo "Omega matter: $OM"
echo "Omega lambda: $OL"
echo "Mstar: $Mstar"
echo "Phistar: $Phistar"
echo "Alfa: $Alfa"
echo "Z low: $zlow"
echo "Z up: $zup"
echo "Z error medio: $zError"
echo "Z desviacion error: $zdError"
echo "Mag lim $magLim"
echo "Delta mag $deltaMag"
echo "Mean error $mError"
echo "Dev error $dError"
echo "Mean Color $meanColor"
echo "Stddev Color $stddevColor"
echo "Mean ErrorColor $meanErrorColor"
echo "Stddev ErrorColor $stddevErrorColor"
echo "área $area"
echo "Salida en $outFile"
echo "estatus: $estatus"
echo "-------------------------"
echo ""

exit($estatus)
