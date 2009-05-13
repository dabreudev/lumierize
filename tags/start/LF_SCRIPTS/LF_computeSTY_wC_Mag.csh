#! /bin/csh

# TITLE: Script para calcular LF con LumFunc_GOYA utilizando STY_wC
# CREATED: 20080320 dabreu
# Pequeño script para calcular la LF utilizando la opción STY_wC en magnitudes
#
# Argumentos [valor por defecto]:
#  - Catálogo de entrada [kk.cat]
#  - Color mean [0.0]
#  - Color stddev [0.0]
#  - magnitud limite [25]
#  - Z up [0]
#  - Z low [0]
#  - área del survey [0.1]
#  - Fichero de salida [kk.lf]
#
# Solo hay dos formas de llamarlo correctamente, dando todos los argumentos o
# ninguno, en cuyo caso tomarán los valores por defecto
#

set estatus=0 #para guardar el status de salida de LumFunc_GOYA

if ($1 == "--help") then
   echo "Uso: LF_computeSTY_wC_Mag.csh inCat colorMean colorStddev magLim zlow zup area outFile"
   echo ""
   echo "Example: LF_computeSTY_wC_Mag.csh kk.cat 0 0.001 25 0 0 0.1 kk.lf"
   exit
endif

if ($1 == "") then
   set inCat="kk.cat"
   set colorMean=0.0
   set colorStddev=0.0
   set magLim=25
   set zlow=0
   set zup=0
   set area=0.1
   set outFile="kk.lf"
else
   set inCat=$1
   set colorMean=$2
   set colorStddev=$3
   set magLim=$4
   set zlow=$5
   set zup=$6
   set area=$7
   set outFile=$8
endif

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
O
75
0.5
f
$inCat
2
7
7
4
4
$colorMean
$colorStddev
m
10
$magLim
25
$zlow
$zup
$area
1
$magLim
$outFile
0
e
COMANDOS

set estatus=$status

echo ""
echo "-------------------------"
echo $0 $*":"
echo "Input file $inCat"
echo "Color mean $colorMean"
echo "Color stddev $colorStddev"
echo "Mag lim $magLim"
echo "Z low: $zlow"
echo "Z up: $zup"
echo "área: $area"
echo "outFile: $outFile"
echo "estatus: $estatus"
echo "-------------------------"
echo ""

exit($estatus)
