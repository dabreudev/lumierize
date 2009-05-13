#! /bin/csh

# TITLE: Script para generar catálogos con LumFunc_GOYA wC
# CREATED: 20080212 dabreu
# Pequeño script para generar un catálogo en magnitudes para realizar
# simulaciones utilizando la opción I de LumFunc_GOYA
#
# Argumentos [valor por defecto]:
#  - Mstar [-20.4]
#  - Phistar [0.0033]
#  - Alfa [-1.3]
#  - Z up [0]
#  - Z low [0]
#  - z error medio [0]
#  - z desviación del error [0]
#  - magnitud limite [25] (salen ~2800 objetos)
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
   echo "Uso: LF_generate_catalog_Mag_wC.csh Mstar Phistar Alfa zlow zup zError zdError magLim mError dError meanColor stddevColor meanErrorColor meanStddevColor area outFile"
   echo ""
   echo "Example: LF_generate_catalog_Mag_wC.csh -20.4 0.0033 -1.3 0 0 0 0 25 0 0 0 0 0 0 0.1 kk.cat"
   exit
endif

if ($1 == "") then
   set Mstar=-20.4
   set Phistar=0.0033
   set Alfa=-1.3
   set zlow=0
   set zup=0
   set zError=0
   set zdError=0
   set magLim=25
   set mError=0.0
   set dError=0.0
   set meanColor=0.0
   set stddevColor=0.0
   set meanErrorColor=0.0
   set stddevErrorColor=0.0
   set area=0.1
   set outFile="kk.cat"
else
   set Mstar=$1
   set Phistar=$2
   set Alfa=$3
   set zlow=$4
   set zup=$5
   set zError=$6
   set zdError=$7
   set magLim=$8
   set mError=$9
   set dError=$10
   set meanColor=$11
   set stddevColor=$12
   set meanErrorColor=$13
   set stddevErrorColor=$14
   set area=$15
   set outFile=$16
endif

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
I
75
0.5
$Mstar
$Phistar
$Alfa
$zlow
$zup
$zError
$zdError
$magLim
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
echo "Mstar: $Mstar"
echo "Phistar: $Phistar"
echo "Alfa: $Alfa"
echo "Z low: $zlow"
echo "Z up: $zup"
echo "Z error medio: $zError"
echo "Z desviacion error: $zdError"
echo "Mag lim $magLim"
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
