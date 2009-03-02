equNo=$1
R1=$2
R2=$3

cat head generatedcode/mid_eq"$equNo"_real tail > generatedcode/vegas_eq"$equNo"_real.cpp
g++ -I. nvegas.o generatedcode/vegas_eq"$equNo"_real.cpp -o exec/vegas_eq"$equNo"_real.exe
exec/vegas_eq"$equNo"_real.exe "$R1" "$R2" > vegasresult/vegas_eq"$equNo"_real.dat

cat head generatedcode/mid_eq"$equNo"_imag tail > generatedcode/vegas_eq"$equNo"_imag.cpp
g++ -I. nvegas.o generatedcode/vegas_eq"$equNo"_imag.cpp -o exec/vegas_eq"$equNo"_imag.exe
exec/vegas_eq"$equNo"_real.exe "$R1" "$R2" > vegasresult/vegas_eq"$equNo"_imag.dat


