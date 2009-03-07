# clean up old file
rm generatedcode/*

# format of the executing: calc_D0Integrand_real.exe <number of equation> p10 p20 p21 p30 p31 p32 m1s m2s m3s m4s
p10=1; p20=5; p21=1; p30=7; p31=15; p32=1
m1s=6561; m2s=8281; m3s=6561; m4s=8281

echo "output integrand equation 1"
exec/calc_D0Integrand_real.exe 1 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq1_real
exec/calc_D0Integrand_imag.exe 1 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq1_imag
echo "output integrand equation 9"
exec/calc_D0Integrand_real.exe 9 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq9_real
exec/calc_D0Integrand_imag.exe 9 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq9_imag

echo "output integrand equation 12"
exec/calc_D0Integrand_real.exe 12 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq12_real
exec/calc_D0Integrand_imag.exe 12 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq12_imag

echo "output integrand equation 18"
# we must pay a little much attention to case equation 18 and on, because it split into z=>0 and z<0
echo "(z>=0)?" > generatedcode/mid_eq18_real
exec/calc_D0Integrand_real.exe 1801 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq18_real
echo ":" >> generatedcode/mid_eq18_real
exec/calc_D0Integrand_real.exe 1802 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq18_real

echo "(z>=0)?" > generatedcode/mid_eq18_imag
exec/calc_D0Integrand_imag.exe 1801 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq18_imag
echo ":" >> generatedcode/mid_eq18_imag
exec/calc_D0Integrand_imag.exe 1802 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq18_imag

echo "output integrand equation 30"
# 4 cases for z, t:
#	z	t	equivalent equation
#	>=0	>=0	30-01
#	>=0	< 0	30-02
#	< 0	>=0	30-03
#	< 0	< 0	30-04

#	REAL CASES
# z>=0, t>=0
echo "(z>=0)?(t>=0?" 	> 	generatedcode/mid_eq30_real
exec/calc_D0Integrand_real.exe 3001 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq30_real
# z>=0, t< 0
echo ":" 		>> 	generatedcode/mid_eq30_real
exec/calc_D0Integrand_real.exe 3002 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq30_real
# z< 0, t>=0
echo "):(t>=0?" 	>> 	generatedcode/mid_eq30_real
exec/calc_D0Integrand_real.exe 3003 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq30_real
# z< 0, t< 0
echo ":" 		>> 	generatedcode/mid_eq30_real
exec/calc_D0Integrand_real.exe 3004 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq30_real
echo ")"		>> 	generatedcode/mid_eq30_real

#	IMAGINARY CASES
# z>=0, t>=0
echo "(z>=0)?(t>=0?" 	> 	generatedcode/mid_eq30_imag
exec/calc_D0Integrand_imag.exe 3001 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq30_imag
# z>=0, t< 0
echo ":" 		>> 	generatedcode/mid_eq30_imag
exec/calc_D0Integrand_imag.exe 3002 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq30_imag
# z< 0, t>=0
echo "):(t>=0?" 	>> 	generatedcode/mid_eq30_imag
exec/calc_D0Integrand_imag.exe 3003 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq30_imag
# z< 0, t< 0
echo ":" 		>> 	generatedcode/mid_eq30_imag
exec/calc_D0Integrand_imag.exe 3004 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq30_imag
echo ")"		>> 	generatedcode/mid_eq30_imag



echo "output integrand equation 39"
# z=>0 and z<0
echo "(z>=0)?" > generatedcode/mid_eq39_real
exec/calc_D0Integrand_real.exe 3901 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq39_real
echo ":" >> generatedcode/mid_eq39_real
exec/calc_D0Integrand_real.exe 3902 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq39_real

echo "(z>=0)?" > generatedcode/mid_eq39_imag
exec/calc_D0Integrand_imag.exe 3901 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq39_imag
echo ":" >> generatedcode/mid_eq39_imag
exec/calc_D0Integrand_imag.exe 3902 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq39_imag




# join the head mid tail to create the code for plotting the function, and the vegas code
# 1 9 12 18 30 39
for equNo in 1 9 12 18 30 39
do	
	echo "building plotting code for integrand ""$equNo"
	cat dataplot4d_head  generatedcode/mid_eq"$equNo"_real  dataplot4d_tail > generatedcode/dataplot_eq"$equNo"_real.cpp
	g++ generatedcode/dataplot_eq"$equNo"_real.cpp -o exec/dataplot_eq"$equNo"_real.exe
	cat dataplot4d_head  generatedcode/mid_eq"$equNo"_imag  dataplot4d_tail > generatedcode/dataplot_eq"$equNo"_imag.cpp
	g++ generatedcode/dataplot_eq"$equNo"_imag.cpp -o exec/dataplot_eq"$equNo"_imag.exe

	echo "building vegas code for integrand ""$equNo"
	cat head generatedcode/mid_eq"$equNo"_real tail > generatedcode/vegas_eq"$equNo"_real.cpp
	g++ -I. nvegas.o generatedcode/vegas_eq"$equNo"_real.cpp -o exec/vegas_eq"$equNo"_real.exe
	cat head generatedcode/mid_eq"$equNo"_imag tail > generatedcode/vegas_eq"$equNo"_imag.cpp
	g++ -I. nvegas.o generatedcode/vegas_eq"$equNo"_imag.cpp -o exec/vegas_eq"$equNo"_imag.exe
done


