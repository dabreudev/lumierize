#! /bin/csh

# TITLE: Script para calcular LF con LumFunc_GOYA utilizando STY
# CREATED: 20060809 dabreu
# Pequeño script para calcular la LF utilizando la opción STY en magnitudes de
# LumFunc_GOYA
#
# MODIFIED: 20060907 dabreu
# Adaptación del cambio de nombres y de parámetros
#
# Argumentos [valor por defecto]:
#  - Catálogo de entrada [kk.cat]
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
   echo "Uso: LF_computeSTY_Mag.csh inCat magLim zlow zup area outFile"
   echo ""
   echo "Example: LF_computeSTY_Mag.csh kk.cat 25 0 0 0.1 kk.lf"
   exit
endif

if ($1 == "") then
   set inCat="kk.cat"
   set magLim=25
   set zlow=0
   set zup=0
   set area=0.1
   set outFile="kk.lf"
else
   set inCat=$1
   set magLim=$2
   set zlow=$3
   set zup=$4
   set area=$5
   set outFile=$6
endif

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
m
75
0.5
f
$inCat
2
7
7
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
echo "Mag lim $magLim"
echo "Z low: $zlow"
echo "Z up: $zup"
echo "área: $area"
echo "outFile: $outFile"
echo "estatus: $estatus"
echo "-------------------------"
echo ""

exit($estatus)
