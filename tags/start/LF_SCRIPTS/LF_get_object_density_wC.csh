#! /bin/csh

# TITLE: Script para obtener densidad de objetos wC con LumFunc_GOYA
# CREATED: 20070214 dabreu
# Pequeño script para obtener la densidad de objetos dados unos parámetros del
# survery utilizando la opción D de LumFunc_GOYA
#
# Argumentos [valor por defecto]:
#  - Magnitudes o Luminosidad [m] (puede ser "m" o "l")
#  - magnitud limite [25]
#  - Mstar [-20.4]
#  - Phistar [0.0033]
#  - Alfa [-1.3]
#  - color_mean [0]
#  - color_stddev [0]
#  - Z up [0]
#  - Z low [0]
#
# Solo hay dos formas de llamarlo correctamente, dando todos los argumentos o
# ninguno, en cuyo caso tomarán los valores por defecto
#

set estatus=0 #para guardar el status de salida de LumFunc_GOYA

if ($1 == "--help") then
   echo "Uso: LF_get_object_density.csh MoL magLim Mstar Phistar Alfa colorMean colorStddev zlow zup"
   exit
endif

if ($1 == "") then
   set MoL="m"
   set magLim=25
   set Mstar=-20.4
   set Phistar=0.0033
   set Alfa=-1.3
   set colorMean=0
   set colorStddev=0.01
   set zlow=0
   set zup=0
else
   set MoL=$1
   set magLim=$2
   set Mstar=$3
   set Phistar=$4
   set Alfa=$5
   set colorMean=$6
   set colorStddev=$7
   set zlow=$8
   set zup=$9
endif

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
d
75
0.5
0
$MoL
$magLim
$Mstar
$Phistar
$Alfa
$colorMean
$colorStddev
$zlow
$zup
e
COMANDOS

set estatus=$status

echo ""
echo "-------------------------"
echo $0 $*":"
echo "Mag o Lum: $MoL"
echo "Mag lim: $magLim"
echo "Mstar: $Mstar"
echo "Phistar: $Phistar"
echo "Alfa: $Alfa"
echo "colorMean: $colorMean"
echo "colorStddev: $colorStddev"
echo "Z low: $zlow"
echo "Z up: $zup"
echo "estatus: $estatus"
echo "-------------------------"
echo ""

exit($estatus)
