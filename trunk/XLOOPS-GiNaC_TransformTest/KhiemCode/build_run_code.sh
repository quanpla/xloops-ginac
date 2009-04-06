p10=1; p20=5; p21=1; p30=7; p31=15; p32=1;
m1s=6561; m2s=8281; m3s=6561; m4s=8281;

#equNo in 1 9 12 18 30 39
equNo=1; dim=4; type=1
#  ./build_run_code.sh $equNo $dim $type $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s

equNo=9; dim=4; type=1
#  ./build_run_code.sh $equNo $dim $type $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s

equNo=12; dim=4; type=1
#  ./build_run_code.sh $equNo $dim $type $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s

equNo=18; dim=3; type=2
#  ./build_run_code.sh $equNo $dim $type $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s

equNo=25; dim=3; type=2
#  ./build_run_code.sh $equNo $dim $type $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s

equNo=30; dim=2; type=4
#  ./build_run_code.sh $equNo $dim $type $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s

equNo=40; dim=1; type=2
./build_run_code.sh $equNo $dim $type $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s

equNo=48; dim=1; type=2
./build_run_code.sh $equNo $dim $type $p10 $p20 $p21 $p30 $p31 $p32 $m1s $m2s $m3s $m4s
