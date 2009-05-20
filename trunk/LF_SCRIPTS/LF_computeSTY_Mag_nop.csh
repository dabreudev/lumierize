#! /bin/csh

# TITLE: Script para calcular LF con LumFunc_GOYA utilizando STY sin poisson
# CREATED: 20090520 dabreu
# Peque�o script para calcular la LF utilizando la opci�n STY en magnitudes de
# LumFunc_GOYA sin la normalizaci�n por poisson
#
# Argumentos [valor por defecto]:
#  - H0 [70]
#  - OM [0.3]
#  - OL [0.7]
#  - Cat�logo de entrada [kk.cat]
#  - magnitud limite [25]
#  - Z up [0]
#  - Z low [0]
#  - �rea del survey [0.1]
#  - Fichero de salida [kk.lf]
#
# Solo hay dos formas de llamarlo correctamente, dando todos los argumentos o
# ninguno, en cuyo caso tomar�n los valores por defecto
#

set estatus=0 #para guardar el status de salida de LumFunc_GOYA

if ($1 == "--help") then
   echo "Uso: LF_computeSTY_Mag.csh H0 OM OL inCat magLim zlow zup area outFile"
   echo ""
   echo "Example: LF_computeSTY_Mag.csh 70 0.3 0.7 kk.cat 25 0 0 0.1 kk.lf"
   exit
endif

if ($1 == "") then
   set H0=70
   set OM=0.3
   set OL=0.7
   set inCat="kk.cat"
   set magLim=25
   set zlow=0
   set zup=0
   set area=0.1
   set outFile="kk.lf"
else
   set H0=$1
   set OM=$2
   set OL=$3
   set inCat=$4
   set magLim=$5
   set zlow=$6
   set zup=$7
   set area=$8
   set outFile=$9
endif

set colz=2
set colm=7

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
m
$H0
$OM
$OL
f
$inCat
$colz
$colm
$colm
m
10
$magLim
25
$zlow
$zup
$area
0
$magLim
$outFile
0
e
COMANDOS

set estatus=$status

echo ""
echo "-------------------------"
echo $0 $*":"
echo "H0 $H0"
echo "Omega matter $OM"
echo "Omega lambda $OL"
echo "Input file $inCat"
echo "Mag lim $magLim"
echo "Z low: $zlow"
echo "Z up: $zup"
echo "�rea: $area"
echo "outFile: $outFile"
echo "estatus: $estatus"
echo "-------------------------"
echo ""

exit($estatus)
