# run the calc_integrand program to create D0 Integrand into mid file, use the same method for image part
# format of the executing: calc_D0Integrand_real.exe <number of equation> p10 p20 p21 p30 p31 p32 m1s m2s m3s m4s
p10=1
p20=5
p21=1
p30=7
p31=15
p32=1
m1s=6561
m2s=8281
m3s=6561
m4s=8281

exec/calc_D0Integrand_real.exe 1 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq1_real
exec/calc_D0Integrand_imag.exe 1 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq1_imag
exec/calc_D0Integrand_imag.exe 9 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq9_imag
exec/calc_D0Integrand_real.exe 9 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq9_real
exec/calc_D0Integrand_real.exe 12 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq12_real
exec/calc_D0Integrand_imag.exe 12 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s > generatedcode/mid_eq12_imag

# we must pay a little much attention to case equation 18, because it split into z=>0 and z<0
echo "(z>=0)?" > generatedcode/mid_eq18_real
exec/calc_D0Integrand_real.exe 1801 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq18_real
echo ":" >> generatedcode/mid_eq18_real
exec/calc_D0Integrand_real.exe 1802 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq18_real

echo "(z>=0)?" > generatedcode/mid_eq18_imag
exec/calc_D0Integrand_imag.exe 1801 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq18_imag
echo ":" >> generatedcode/mid_eq18_imag
exec/calc_D0Integrand_imag.exe 1802 $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s >> generatedcode/mid_eq18_imag


# join the head mid tail to create the code for plotting the function
# 	input for sub_generate_plotdata.sh : <no of D0 Integrand equation> <start> <step> <end>
# see the result in plotdata folder
start=-1000
step=1
end=1000
./sub_generate_plotdata.sh 1 $start $step $end
./sub_generate_plotdata.sh 9 $start $step $end
./sub_generate_plotdata.sh 12 $start $step $end
./sub_generate_plotdata.sh 18 $start $step $end

# join the head mid tail to create the code for vegas
# see the result in vegasresult folder
R1 = -1000
R2 = 1000
./sub_run_vegas.sh 1 $R1 $R2
./sub_run_vegas.sh 9 $R1 $R2
./sub_run_vegas.sh 12 $R1 $R2
./sub_run_vegas.sh 18 $R1 $R2

