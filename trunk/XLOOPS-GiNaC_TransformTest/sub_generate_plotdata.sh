equNo=$1
first=$2
step=$3
last=$4

cat dataplot4d_head generatedcode/mid_eq"$equNo"_real dataplot4d_tail > generatedcode/dataplot_eq"$equNo"_real.cpp
g++ generatedcode/dataplot_eq"$equNo"_real.cpp -o exec/dataplot_eq"$equNo"_real.exe
exec/dataplot_eq"$equNo"_real.exe "$first" "$step" "$last" > plotdata/eq"$equNo"_real.dat

cat dataplot4d_head generatedcode/mid_eq"$equNo"_imag dataplot4d_tail > generatedcode/dataplot_eq"$equNo"_imag.cpp
g++ generatedcode/dataplot_eq"$equNo"_imag.cpp -o exec/dataplot_eq"$equNo"_imag.exe
exec/dataplot_eq"$equNo"_imag.exe "$first" "$step" "$last" > plotdata/eq"$equNo"_imag.dat

