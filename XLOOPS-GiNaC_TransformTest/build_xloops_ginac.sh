# Build underly XLOOPS-GiNaC programs
g++ `pkg-config --cflags --libs ginac` -c my_fns.cpp
g++ `pkg-config --cflags --libs ginac` -c lev1.cpp
g++ `pkg-config --cflags --libs ginac` -c D0Integrand.cpp

# Build the executable
g++ `pkg-config --cflags --libs ginac` my_fns.o lev1.o D0Integrand.o calc_D0Integrand_real.cpp -o calc_D0Integrand_real.exe
g++ `pkg-config --cflags --libs ginac` my_fns.o lev1.o D0Integrand.o calc_D0Integrand_imag.cpp -o calc_D0Integrand_imag.exe

