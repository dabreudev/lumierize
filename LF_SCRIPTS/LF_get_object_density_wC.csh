#! /bin/csh

# TITLE: Script para obtener densidad de objetos wC con LumFunc_GOYA
# CREATED: 20070214 dabreu
# Peque�o script para obtener la densidad de objetos dados unos par�metros del
# survery utilizando la opci�n D de LumFunc_GOYA
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
#  - color_mean [0]
#  - color_stddev [0]
#  - Z up [0]
#  - Z low [0]
#
# Solo hay dos formas de llamarlo correctamente, dando todos los argumentos o
# ninguno, en cuyo caso tomar�n los valores por defecto
#

set estatus=0 #para guardar el status de salida de LumFunc_GOYA

if ($1 == "--help") then
   echo "Uso: LF_get_object_density.csh H0 OM OL MoL magLim Mstar Phistar Alfa colorMean colorStddev zlow zup"
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
   set colorMean=0
   set colorStddev=0.01
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
   set colorMean=$9
   set colorStddev=$10
   set zlow=$11
   set zup=$12
endif

$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
d
$H0
$OM
$OL
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
echo "H0: $H0"
echo "OM: $OM"
echo "OL: $OL"
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
