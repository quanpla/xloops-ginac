#	1.	Get input range
if test $1
then
	# input avaialable
	R1=$1; R2=$2;
else 
	# edit this to change the default limit range of the integral
	R1=-3000; R2=3000; # -3000..3000 will submerge the place with most anomalies in D0integral -> best practice
fi

#	2.	Sset default output folder
outputFolder="vegasresult/"

# 	3.	The main loop (loop all equation)
for equNo in 75
do
	# generate real  part
	outputFileName="NPoint_vegas_equ""$equNo""_real""_$R1""_$R2"".dat"
	echo "Calculating vegas integral for integrand $equNo - real  part. File name: ""$outputFolder""$outputFileName"
	exec/NPoint_vegas_eq"$equNo"_real.exe "$R1" "$R2" > "$outputFolder""$outputFileName"

	# generate image part
	outputFileName="NPoint_vegas_equ""$equNo""_imag""_$R1""_$R2"".dat"
	echo "Calculating vegas integral for integrand $equNo - imag  part. File name: ""$outputFolder""$outputFileName"
	exec/NPoint_vegas_eq"$equNo"_imag.exe "$R1" "$R2" > "$outputFolder""$outputFileName"

#end of looping for each integrand equation
done

