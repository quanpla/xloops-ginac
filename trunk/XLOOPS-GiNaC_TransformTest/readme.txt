
Step #		Command				Description
1		./build_xloops_ginac.sh		run this to rebuild the XLOOPS-GiNaC code

2		./build_run_code.sh		run this to build all generated code (result in exec folder)

3		./run_plotdata.sh		Generate data to plot
						Result in folder plotdata/
						It has parameter: first step last. If parameter not provided, it will use 							default	values 

4		./run_vegas.sh			Executing vegas Integral calculation
						Result in folder vegasresult
						It has parameter: R1 R2 (the integral limit). If parameter not provided, it 							will use default values 


for a quick reference:
all generated code will go to generatedcode folder
all executable code will go to exec folder
all data file (to plot) should go to plotdata folder
all vegas executing result should go to vegasresult folder
