#! /bin/csh

# TITLE: Script para calcular LF con LumFunc_GOYA utilizando STY_gmz_M_wC
# CREATED: 20090319 dabreu
# Pequeño script para calcular la LF utilizando la opción STY_gmz_M_wC
#
# Argumentos [valor por defecto]:
#  - Catálogo de entrada [kk.cat]
#  - Color mean [0.0]
#  - Color stddev [0.0]
#  - magnitud limite [25]
#  - deltamag de funcion de seleccion [0.5]
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
   echo "Uso: LF_computeSTY_gmz_M_wC.csh inCat colorMean colorStddev magLim deltmag_fsel zlow zup area outFile"
   echo ""
   echo "Example: LF_computeSTY_gmz_M_wC.csh kk.cat 0 0.001 25 0 0 0 0.1 kk.lf"
   exit
endif

if ($1 == "") then
   set inCat="kk.cat"
   set colorMean=0.0
   set colorStddev=0.0
   set magLim=25
   set deltamag=0.5
   set zlow=0
   set zup=0
   set area=0.1
   set outFile="kk.lf"
else
   set inCat=$1
   set colorMean=$2
   set colorStddev=$3
   set magLim=$4
   set deltamag=$5
   set zlow=$6
   set zup=$7
   set area=$8
   set outFile=$9
endif

set colz=2
set colerrz=3
set colmSel=7
set colmDist=12
set colerrmDist=8

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
R
75
0.5
f
$inCat
$colz
$colerrz
$colmSel
$colmDist
$colerrmDist
$colorMean
$colorStddev
m
10
$magLim
25
$zlow
$zup
$area
$magLim
$deltamag
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
echo "Delta mag (Fermi) $deltamag"
echo "Z low: $zlow"
echo "Z up: $zup"
echo "área: $area"
echo "outFile: $outFile"
echo "estatus: $estatus"
echo "-------------------------"
echo ""

exit($estatus)
