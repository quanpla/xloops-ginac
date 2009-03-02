# Build underly XLOOPS-GiNaC programs
g++ `pkg-config --cflags --libs ginac` -c my_fns.cpp
g++ `pkg-config --cflags --libs ginac` -c trmchk.cpp
g++ `pkg-config --cflags --libs ginac` -c trm2F.cpp
g++ `pkg-config --cflags --libs ginac` -c lev1.cpp
g++ `pkg-config --cflags --libs ginac` -c lev2.cpp
g++ `pkg-config --cflags --libs ginac` -c lev3.cpp
g++ `pkg-config --cflags --libs ginac` -c D0Integrand.cpp

# Build the executable and move to exec folder
g++ `pkg-config --cflags --libs ginac` my_fns.o lev1.o lev2.o lev3.o trm2F.o trmchk.o D0Integrand.o calc_D0Integrand_real.cpp -o exec/calc_D0Integrand_real.exe
g++ `pkg-config --cflags --libs ginac` my_fns.o lev1.o lev2.o lev3.o trm2F.o trmchk.o D0Integrand.o calc_D0Integrand_imag.cpp -o exec/calc_D0Integrand_imag.exe

# build vegas before hand
gcc -c nvegas.c

