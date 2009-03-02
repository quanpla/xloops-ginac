# run the calc_integrand program to create D0 Integrand into mid file, use the same method for image part
# format of the executing: calc_D0Integrand_real.exe <number of equation> p10 p20 p21 p30 p31 p32 m1s m2s m3s m4s
./calc_D0Integrand_real.exe 1 1 5 1 7 15 1 6561 8281 6561 8281 > generatedcode/mid_eq1_real
./calc_D0Integrand_real.exe 9 1 5 1 7 15 1 6561 8281 6561 8281 > generatedcode/mid_eq9_real
./calc_D0Integrand_real.exe 12 1 5 1 7 15 1 6561 8281 6561 8281 > generatedcode/mid_eq12_real
./calc_D0Integrand_real.exe 18 1 5 1 7 15 1 6561 8281 6561 8281 > generatedcode/mid_eq18_real

# join the head mid tail to create the code for plotting the function
cat dataplot4d_head generatedcode/mid_eq1_real dataplot4d_tail > generatedcode/dataplot_eq1_real.cpp
# please try to make follwing script run
g++ generatedcode/dataplot_eq1_real.cpp -o generatedcode/dataplot_eq1_real.exe
generatedcode/dataplot_eq1_real.exe > plotdata/eq1_real.dat

# as well as making following script run
# join the head mid tail to create the code for vegas
cat head mid_eq1_real tail > generatedcode/vegas_eq1_real.cpp
g++ nvegas.o generatedcode/vegas_eq1_real.cpp -o generatedcode/vegas_eq1_real.exe
