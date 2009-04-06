# Build underly XLOOPS-GiNaC programs

for fileName in "my_fns" "trmchk" "trm2F" "lev1" "lev2" "lev3" "lev4" "D0Integrand"
do
	if ! [ -f "$fileName.o" ];
	then
		echo "Building the object file $fileName.o"
		g++ `pkg-config --cflags --libs ginac` -c "$fileName".cpp
	else
			echo "$fileName.o existed! No thing happened."
	fi
done

# Build the executable and move to exec folder
echo "building exec/calc_D0Integrand_real.exe"
g++ `pkg-config --cflags --libs ginac` my_fns.o lev1.o lev2.o lev3.o lev4.o trm2F.o trmchk.o D0Integrand.o calc_D0Integrand_real.cpp -o exec/calc_D0Integrand_real.exe

echo "building exec/calc_D0Integrand_imag.exe"
g++ `pkg-config --cflags --libs ginac` my_fns.o lev1.o lev2.o lev3.o lev4.o trm2F.o trmchk.o D0Integrand.o calc_D0Integrand_imag.cpp -o exec/calc_D0Integrand_imag.exe

# build vegas before hand
if ! [ -f "nvegas.o" ];
then
	echo "building nvegas.o"
	gcc -c nvegas.c
else
	echo "nvegas.o existed! No thing happened."
fi
