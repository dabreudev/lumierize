#! /bin/csh

# TITLE: Script para obtener densidad de objetos con LumFunc_GOYA
# CREATED: 20060905 dabreu
# Pequeño script para obtener la densidad de objetos dados unos parámetros del
# survery utilizando la opción C de LumFunc_GOYA
#
# MODIFIED: 20090518 dabreu
# Support for new cosmology
#
# Argumentos [valor por defecto]:
#  - H0 [70]
#  - Omega matter [0.3]
#  - Omega lambda [0.7]
#  - Magnitudes o Luminosidad [m] (puede ser "m" o "l")
#  - magnitud limite [25]
#  - Mstar [-20.4]
#  - Phistar [0.0033]
#  - Alfa [-1.3]
#  - Z up [0]
#  - Z low [0]
#
# Solo hay dos formas de llamarlo correctamente, dando todos los argumentos o
# ninguno, en cuyo caso tomarán los valores por defecto
#

set estatus=0 #para guardar el status de salida de LumFunc_GOYA

if ($1 == "--help") then
   echo "Uso: LF_get_object_density.csh H0 OM OL MoL magLim Mstar Phistar Alfa zlow zup"
   echo ""
   echo "Example: LF_get_object_density.csh 70 0.3 0.7 m 25 -20.4 0.0033 -1.3 0 0"
   exit
endif

if ($1 == "") then
   set H0=70
   set OM=0.3
   set OL=0.7
   set MoL="m"
   set magLim=25
   set Mstar=-20.4
   set Phistar=0.0033
   set Alfa=-1.3
   set zlow=0
   set zup=0
else
   set H0=$1
   set OM=$2
   set OL=$3
   set MoL=$4
   set magLim=$5
   set Mstar=$6
   set Phistar=$7
   set Alfa=$8
   set zlow=$9
   set zup=$10
endif

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
c
$H0
$OM
$OL
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
echo "H0: $H0"
echo "Omega matter: $OM"
echo "Omega lambda: $OL"
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
