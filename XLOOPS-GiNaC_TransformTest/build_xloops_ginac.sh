# Build underly XLOOPS-GiNaC programs
echo "building my_fns.o"
g++ `pkg-config --cflags --libs ginac` -c my_fns.cpp

echo "building trmchk.o"
g++ `pkg-config --cflags --libs ginac` -c trmchk.cpp

echo "building trm2F.o"
g++ `pkg-config --cflags --libs ginac` -c trm2F.cpp

echo "building lev1.o"
g++ `pkg-config --cflags --libs ginac` -c lev1.cpp

echo "building lev2.o"
g++ `pkg-config --cflags --libs ginac` -c lev2.cpp

echo "building lev3.o"
g++ `pkg-config --cflags --libs ginac` -c lev3.cpp

echo "building D0Integrand.o"
g++ `pkg-config --cflags --libs ginac` -c D0Integrand.cpp

# Build the executable and move to exec folder
echo "building exec/calc_D0Integrand_real.exe"
g++ `pkg-config --cflags --libs ginac` my_fns.o lev1.o lev2.o lev3.o trm2F.o trmchk.o D0Integrand.o calc_D0Integrand_real.cpp -o exec/calc_D0Integrand_real.exe

echo "building exec/calc_D0Integrand_imag.exe"
g++ `pkg-config --cflags --libs ginac` my_fns.o lev1.o lev2.o lev3.o trm2F.o trmchk.o D0Integrand.o calc_D0Integrand_imag.cpp -o exec/calc_D0Integrand_imag.exe

# build vegas before hand
echo "building nvegas.o"
gcc -c nvegas.c

