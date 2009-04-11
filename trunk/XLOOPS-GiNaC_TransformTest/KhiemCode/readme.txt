I.	Execution Plan & Description

============================================================================================================
Step #		Command				Description
========	========================	============================================================
1		./build_xloops_ginac.sh		run this to rebuild the XLOOPS-GiNaC code
					
2		./build_run_code.sh		run this to build all generated code (result in exec folder)
					
3		./run_plotdata.sh		Generate data to plot
						Result in folder plotdata/
						It has parameter: first step last.
						If parameter not provided, it will use default values 
					
4		./run_vegas.sh			Executing vegas Integral calculation
						Result in folder vegasresult
						It has parameter: R1 R2 (the integral limit).
						If parameter not provided, it will use default values 
============================================================================================================
		Table	1. Description for script files.


II.	Sample Execution

$>./build_xloops_ginac.sh
$>./build_run_code.sh
$>./run_plotdata.sh -1000 1 1000
$>./run_vegas.sh -1000 1000


III.	Sample Script File for Batch Run of Many Integral Limits:

$>vim batch_run_vegas.sh
for limit in 100 1000 10000 100000 1000000
do
	./run_vegas.sh -$limit $limit
done
$>chmod +x batch_run_vegas.sh
$>./batch_run_vegas.sh


IV.	Folder Description

For a quick reference
all generated code will go to generatedcode folder
all executable code will go to exec folder
all data file (to plot) should go to plotdata folder
all vegas executing result should go to vegasresult folder

IV.	Howto add a new build_run_code
The script build_run_code calls build_run_1code to create:
1./ The source file
	+
2./ The exe file
