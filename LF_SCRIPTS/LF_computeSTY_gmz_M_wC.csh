#! /bin/csh

# TITLE: Script para calcular LF con LumFunc_GOYA utilizando STY_gmz_M_wC
# CREATED: 20090319 dabreu
# Pequeño script para calcular la LF utilizando la opción STY_gmz_M_wC
#
# MODIFIED: 20090518 dabreu
# Support for new cosmology
#
# Argumentos [valor por defecto]:
#  - H0 [70]
#  - Omega matter [0.3]
#  - Omega lambda [0.7]
#  - Catálogo de entrada [kk.cat]
#  - Color mean [0.0]
#  - Color stddev [0.0]
#  - magnitud limite [25]
#  - deltamag de funcion de seleccion [0.0]
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
   echo "Uso: LF_computeSTY_gmz_M_wC.csh H0 OM OL inCat colorMean colorStddev magLim deltmag_fsel zlow zup area outFile"
   echo ""
   echo "Example: LF_computeSTY_gmz_M_wC.csh 70 0.3 0.7 kk.cat 0 0.001 25 0.0 0 0 0.1 kk.lf"
   exit
endif

if ($1 == "") then
   set H0=70
   set OM=0.3
   set OL=0.7
   set inCat="kk.cat"
   set colorMean=0.0
   set colorStddev=0.0
   set magLim=25
   set deltamag=0.0
   set zlow=0
   set zup=0
   set area=0.1
   set outFile="kk.lf"
else
   set H0=$1
   set OM=$2
   set OL=$3
   set inCat=$4
   set colorMean=$5
   set colorStddev=$6
   set magLim=$7
   set deltamag=$8
   set zlow=$9
   set zup=$10
   set area=$11
   set outFile=$12
endif

set colz=2
set colerrz=3
set colmSel=7
set colmDist=12
set colerrmDist=13

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
R
$H0
$OM
$OL
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
echo "H0 $H0"
echo "Omega matter $OM"
echo "Omega lambda $OL"
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
