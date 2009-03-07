#	1.	Get input range
if test $1
then
	# input avaialable
	firstNum=$1; step=$2; lastNum=$3;
else 
	# edit this to change the default range of plotting data
	firstNum=-100; step=1; lastNum=100;
fi

#	2.	Sset default output folder
outputFolder="plotdata/";

# 	3.	The main loop (loop all dataplot)
for equNo in 1 9 12 18 30 39
do
	# generate real  part
	outputFileName="eq""$equNo""_real""_$firstNum""_$step""_$lastNum"".dat";
	echo "Generating plotting data for integrand $equNo - real  part. File name: ""$outputFolder""$outputFileName";
	exec/dataplot_eq"$equNo"_real.exe $firstNum $step $lastNum > "$outputFolder""$outputFileName"

	# generate image part
	outputFileName="eq""$equNo""_imag""_$firstNum""_$step""_$lastNum"".dat";
	echo "Generating plotting data for integrand $equNo - imag  part. File name: ""$outputFolder""$outputFileName";
	exec/dataplot_eq"$equNo"_imag.exe $firstNum $step $lastNum > "$outputFolder""$outputFileName";

#end of looping for each integrand equation
done

