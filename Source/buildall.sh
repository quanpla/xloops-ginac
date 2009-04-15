for fileName in "my_fns" "trmchk" "lev1" "lev2" "lev3" "lev4" "lev5"
do
        if ! [ -f "$fileName.o" ];
        then
                echo "Building the object file $fileName.o"
                g++ `pkg-config --cflags --libs ginac` -c "$fileName".cpp
        else
                        echo "$fileName.o existed! No thing happened."
        fi
done

echo "building test.exe"
g++ `pkg-config --cflags --libs ginac` -o test.exe test.cpp my_fns.o trmchk.o lev1.o lev2.o lev3.o lev4.o lev5.o

echo "now execute: ./test.exe OR ./run.sh"

