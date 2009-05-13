#! /bin/csh

# TITLE: Script para obtener densidad de objetos con LumFunc_GOYA
# CREATED: 20060905 dabreu
# Peque�o script para obtener la densidad de objetos dados unos par�metros del
# survery utilizando la opci�n C de LumFunc_GOYA
#
# Argumentos [valor por defecto]:
#  - Magnitudes o Luminosidad [m] (puede ser "m" o "l")
#  - magnitud limite [25]
#  - Mstar [-20.4]
#  - Phistar [0.0033]
#  - Alfa [-1.3]
#  - Z up [0]
#  - Z low [0]
#
# Solo hay dos formas de llamarlo correctamente, dando todos los argumentos o
# ninguno, en cuyo caso tomar�n los valores por defecto
#

set estatus=0 #para guardar el status de salida de LumFunc_GOYA

if ($1 == "--help") then
   echo "Uso: LF_get_object_density.csh MoL magLim Mstar Phistar Alfa zlow zup"
   echo ""
   echo "Example: LF_get_object_density.csh m 25 -20.4 0.0033 -1.3 0 0"
   exit
endif

if ($1 == "") then
   set MoL="m"
   set magLim=25
   set Mstar=-20.4
   set Phistar=0.0033
   set Alfa=-1.3
   set zlow=0
   set zup=0
else
   set MoL=$1
   set magLim=$2
   set Mstar=$3
   set Phistar=$4
   set Alfa=$5
   set zlow=$6
   set zup=$7
endif

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
c
75
0.5
0
$MoL
$magLim
$Mstar
$Phistar
$Alfa
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
echo "Z low: $zlow"
echo "Z up: $zup"
echo "estatus: $estatus"
echo "-------------------------"
echo ""

exit($estatus)
