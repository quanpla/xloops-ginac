1./ Quick notes:
To edit what to print to the screen, edit test.cpp, for example:
	cout << fn_alpha(l, k);
For D0, no need for include "Lev1", Lev2..., see the file testD0.cpp

2./ File description:

File            Description
___________     ___________________________________________________
buildall.sh     build all xloops-ginac file
                build test.exe

cleanup.sh      clean *.o and *~ files

extrm.h         define the variables for extern
                (for programming purpose, you merely edit)

lev1->lev5      functions to calculate intermidiate terms
                I have changed the way to calculate
                Have a look at it

my_fns          functions such as my_step, my_csgn, etc.

trmchk          check the validity, e.g. denom <> 0

trm.h           define the global variables

test.cpp        The file that contains all PRINTABLE info
                edit it as your pleasure.
testD0.cpp	Only one call for D0, easy do it

So you can edit: buildall.sh, test.cpp, testD0.cpp
