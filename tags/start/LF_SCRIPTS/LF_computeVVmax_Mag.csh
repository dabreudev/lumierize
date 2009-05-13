#! /bin/csh

# TITLE: Script para calcular LF con LumFunc_GOYA utilizando VVMAX
# CREATED: 20060425 dabreu
# Pequeño script para calcular la LF utilizando VVMAX de LumFunc_GOYA en
# Magnitudes.
#
# MODIFIED: 20060905 dabreu
# Ajustes con los cambios de nombres y algunos parámetros
# MODIFIED: 20061212 dabreu
# Ajustes con los cambios de algunos parámetros
# MODIFIED: 20080317 dabreu
# New param: File for full output.
#
# Argumentos [valor por defecto]:
#  - Catálogo de entrada [kk.cat]
#  - Z up [0]
#  - Z low [0]
#  - área del survey [0.1]
#  - magnitud limite [25]
#  - nombre fichero de salida [kk.lf]
#  - Full output fileName [kk_full.lf]
#
# Solo hay dos formas de llamarlo correctamente, dando todos los argumentos o
# ninguno, en cuyo caso tomarán los valores por defecto
#

set estatus=0 #para guardar el status de salida de LumFunc_GOYA

if ($1 == "--help") then
   echo "Uso: LF_computeVVmax_Mag.csh inCat zlow zup area magLim outFile outFullFile"
   echo ""
   echo "Example: LF_computeVVmax_Mag.csh kk.cat 0 0 0.1 25 kk.lf kk_full.lf"
   exit
endif

if ($1 == "") then
   set inCat="kk.cat"
   set zlow=0
   set zup=0
   set area=0.1
   set magLim=25
   set outFile="kk.lf"
   set outFullFile = "kk_full.lf"
else
   set inCat=$1
   set zlow=$2
   set zup=$3
   set area=$4
   set magLim=$5
   set outFile=$6
   set outFullFile=$7
endif

# Hay que poner el link a estable cuando LumFunc_GOYA no esté ocupado
# y cambiar las columnas para que utilice los colores (no será siempre 7 7 7 7)
$OPERA_INST/bin/LumFunc_GOYA << COMANDOS
v
75
0.5
f
$inCat
2
7
7
7
7
m
-30
-8
25
$zlow
$zup
$area
0
$magLim
$outFile
1
$outFullFile
e
COMANDOS

set estatus=$status

echo ""
echo "-------------------------"
echo $0 $*":"
echo "Input file $inCat"
echo "Z low $zlow"
echo "Z up $zup"
echo "area $area"
echo "Mag lim $magLim"
echo "Out file $outFile"
echo "Out full file $outFullFile"
echo "estatus: $estatus"
echo "-------------------------"
echo ""

exit($estatus)
