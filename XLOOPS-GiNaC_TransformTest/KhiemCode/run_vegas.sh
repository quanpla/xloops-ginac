#	1.	Get input range
if test $1
then
	# input avaialable
	R1=$1; R2=$2;
else 
	# edit this to change the default limit range of the integral
	R1=-100000; R2=100000;
fi

#	2.	Sset default output folder
outputFolder="vegasresult/"

# 	3.	The main loop (loop all equation)
#1 9 12 18 25 30 40 48
for equNo in 40 48
do
	# generate real  part
	outputFileName="vegas_equ""$equNo""_real""_$R1""_$R2"".dat"
	echo "Calculating vegas integral for integrand $equNo - real  part. File name: ""$outputFolder""$outputFileName"
	exec/vegas_eq"$equNo"_real.exe "$R1" "$R2" > "$outputFolder""$outputFileName"

	# generate image part
	outputFileName="vegas_equ""$equNo""_imag""_$R1""_$R2"".dat"
	echo "Calculating vegas integral for integrand $equNo - imag  part. File name: ""$outputFolder""$outputFileName"
	exec/vegas_eq"$equNo"_imag.exe "$R1" "$R2" > "$outputFolder""$outputFileName"

#end of looping for each integrand equation
done

