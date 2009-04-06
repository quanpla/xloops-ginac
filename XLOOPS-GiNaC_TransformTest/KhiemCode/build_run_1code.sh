# get inputs
equNo=$1;
dim=$2;
type=$3;
# type=1 means trivial case integrand like integrand 1.
# type=2 means 2 cases like integrand 18
# type=4 means 4 cases like integrand 30

shift 3;
p10=$1; p20=$2; p21=$3; p30=$4; p31=$5; p32=$6;

shift 6;
m1s=$1; m2s=$2; m3s=$3; m4s=$4;

echo $equNo $dim $type $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s

echo "output integrand equation $equNo";

if [ $type -eq 1 ];
then
	exec/calc_D0Integrand_real.exe $equNo $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq"$equNo"_real;
	exec/calc_D0Integrand_imag.exe $equNo $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq"$equNo"_imag;
elif [ $type -eq 2 ];
then
# we must pay a little much attention to case equation 18 and on, because it split into z=>0 and z<0
	echo "(z>=0)?" > generatedcode/mid_eq"$equNo"_real;
	exec/calc_D0Integrand_real.exe "$equNo"01 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq"$equNo"_real;
	echo ":" >> generatedcode/mid_eq"$equNo"_real;
	exec/calc_D0Integrand_real.exe "$equNo"02 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq"$equNo"_real;

	echo "(z>=0)?" > generatedcode/mid_eq"$equNo"_imag;
	exec/calc_D0Integrand_imag.exe "$equNo"01 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq"$equNo"_imag;
	echo ":" >> generatedcode/mid_eq"$equNo"_imag;
	exec/calc_D0Integrand_imag.exe "$equNo"02 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq"$equNo"_imag;
elif [ $type -eq 4 ];
then
#	z	t	equivalent equation
#	>=0	>=0	$equNo-01
#	>=0	< 0	$equNo-02
#	< 0	>=0	$equNo-03
#	< 0	< 0	$equNo-04

#	REAL CASES
# z>=0, t>=0
	echo "(z>=0)?(t>=0?" 	> 	generatedcode/mid_eq$equNo_real;
	exec/calc_D0Integrand_real.exe "$equNo"01 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq$equNo_real;
# z>=0, t< 0
	echo ":" 		>> 	generatedcode/mid_eq$equNo_real;
	exec/calc_D0Integrand_real.exe "$equNo"02 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq$equNo_real;
# z< 0, t>=0
	echo "):(t>=0?" 	>> 	generatedcode/mid_eq$equNo_real;
	exec/calc_D0Integrand_real.exe "$equN"o03 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq$equNo_real;
# z< 0, t< 0
	echo ":" 		>> 	generatedcode/mid_eq$equNo_real;
	exec/calc_D0Integrand_real.exe "$equNo"04 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq$equNo_real;
	echo ")"		>> 	generatedcode/mid_eq$equNo_real;

#	IMAGINARY CASES
# z>=0, t>=0
	echo "(z>=0)?(t>=0?" 	> 	generatedcode/mid_eq$equNo_imag
	exec/calc_D0Integrand_imag.exe "$equNo"01 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq$equNo_imag;
# z>=0, t< 0
	echo ":" 		>> 	generatedcode/mid_eq$equNo_imag;
	exec/calc_D0Integrand_imag.exe "$equNo"02 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq$equNo_imag;
# z< 0, t>=0
	echo "):(t>=0?" 	>> 	generatedcode/mid_eq$equNo_imag;
	exec/calc_D0Integrand_imag.exe "$equNo"03 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq$equNo_imag;
# z< 0, t< 0
	echo ":" 		>> 	generatedcode/mid_eq$equNo_imag;
	exec/calc_D0Integrand_imag.exe "$equNo"04 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq$equNo_imag;
	echo ")"		>> 	generatedcode/mid_eq$equNo_imag;
fi



# join the head mid tail to create the code for plotting the function, and the vegas code

echo "building plotting code for integrand $equNo";
cat headtail/dataplot"$dim"d_head  generatedcode/mid_eq"$equNo"_real  headtail/dataplot"$dim"d_tail > generatedcode/dataplot_eq"$equNo"_real.cpp;
g++ generatedcode/dataplot_eq"$equNo"_real.cpp -o exec/dataplot_eq"$equNo"_real.exe
cat headtail/dataplot"$dim"d_head  generatedcode/mid_eq"$equNo"_imag  headtail/dataplot"$dim"d_tail > generatedcode/dataplot_eq"$equNo"_imag.cpp
g++ generatedcode/dataplot_eq"$equNo"_imag.cpp -o exec/dataplot_eq"$equNo"_imag.exe

echo "building vegas code for integrand $equNo"
cat headtail/head"$dim" generatedcode/mid_eq"$equNo"_real headtail/tail > generatedcode/vegas_eq"$equNo"_real.cpp
g++ -I. ./nvegas.o generatedcode/vegas_eq"$equNo"_real.cpp -o exec/vegas_eq"$equNo"_real.exe
cat headtail/head"$dim" generatedcode/mid_eq"$equNo"_imag headtail/tail > generatedcode/vegas_eq"$equNo"_imag.cpp
g++ -I. ./nvegas.o generatedcode/vegas_eq"$equNo"_imag.cpp -o exec/vegas_eq"$equNo"_imag.exe
