#!/bin/csh

echo Enter directory where the tar was untared:
#setenv untared_dir $<
setenv untared_dir /home/dabreu/LF/OPERA
setenv install_dir $untared_dir/
setenv modulos_dir $untared_dir/modulos
setenv include_dir $untared_dir/include
echo "Enter directory to install binaries (this must be in your PATH)"
#setenv bin_dir $<
setenv bin_dir $install_dir/bin
echo $bin_dir > variable_bin_dir
echo "Enter directory to install libraries"
#setenv lib_dir $<
setenv lib_dir $install_dir/lib
echo $lib_dir > variable_lib_dir
echo $untared_dir > variable_untared_dir
echo PASO
if($?PGPLOT_DIR) then
else
  echo PGPLOT_DIR variable is not set. It must point to PGPLOT   Library:
  setenv PGPLOT_DIR $<
endif
if($?FITSIO_DIR) then
else
  echo FITSIO_DIR variable is not set. It must point to FITSIO   Library:
  setenv FITSIO_DIR $<
endif
if($?WCSTOOLS_DIR) then
else
  echo WCSTOOLS_DIR variable is not set. It must point to WCSTOOLS Library:
  setenv WCSTOOLS_DIR $<
endif
if($?CBUTTON_DIR) then
else
  echo CBUTTON_DIR variable is not set. It must point to C Button Library:
  setenv CBUTTON_DIR $<
endif

echo 'Which operating system are you using?'
echo    1   Linux
echo    2   RedHat 6.x
echo    3   RedHat 7.x
echo    4   Digital Unix
echo    5   SunOS
echo    6   Solaris
setenv OS $<

switch ($OS)

case 1:
        echo Compiling for Linux
	cd $untared_dir
	echo "OPERA_DIR = $untared_dir" > Makefile
	echo "BIN_OPERA_DIR = $bin_dir" >> Makefile
	echo "LIB_OPERA_DIR = $lib_dir" >> Makefile
	cat $untared_dir/Makefile.sample >> Makefile
	make linux
 	breaksw
case 2:
        echo Compiling for RedHat 6.x
	cd $untared_dir
	echo "OPERA_DIR = $untared_dir" > Makefile
	echo "BIN_OPERA_DIR = $bin_dir" >> Makefile
	echo "LIB_OPERA_DIR = $lib_dir" >> Makefile
	cat $untared_dir/Makefile.sample >> Makefile
	make redhat6
 	breaksw
case 3:
        echo Compiling for RedHat 7.x
	cd $untared_dir
	echo "OPERA_DIR = $untared_dir" > Makefile
	echo "BIN_OPERA_DIR = $bin_dir" >> Makefile
	echo "LIB_OPERA_DIR = $lib_dir" >> Makefile
	cat $untared_dir/Makefile.sample >> Makefile
	make redhat7
 	breaksw

case 4:
        echo Compiling for Digital Unix
	cd $untared_dir
	echo "OPERA_DIR = $untared_dir" > Makefile
	echo "BIN_OPERA_DIR = $bin_dir" >> Makefile
	echo "LIB_OPERA_DIR = $lib_dir" >> Makefile
	cat $untared_dir/Makefile.sample >> Makefile
	make digital
 	breaksw
case 5:
        echo Compiling for SunOS
	cd $untared_dir
	echo "OPERA_DIR = $untared_dir" > Makefile
	echo "BIN_OPERA_DIR = $bin_dir" >> Makefile
	echo "LIB_OPERA_DIR = $lib_dir" >> Makefile
	cat $untared_dir/Makefile.sample >> Makefile
	make sun
 	breaksw
case 6:
        echo Compiling for Solaris
	cd $untared_dir
	echo "OPERA_DIR = $untared_dir" > Makefile
	echo "BIN_OPERA_DIR = $bin_dir" >> Makefile
	echo "LIB_OPERA_DIR = $lib_dir" >> Makefile
	cat $untared_dir/Makefile.sample >> Makefile
	make solaris
 	breaksw
endsw
exit


echo Installation finished
echo OPERA ready to use



