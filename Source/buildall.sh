g++ -c `pkg-config --cflags --libs ginac` trm2F.cpp
g++ -c `pkg-config --cflags --libs ginac` my_fns.cpp
g++ -c `pkg-config --cflags --libs ginac` trmchk.cpp
g++ -c `pkg-config --cflags --libs ginac` lev1.cpp
g++ -c `pkg-config --cflags --libs ginac` lev2.cpp
g++ -c `pkg-config --cflags --libs ginac` lev3.cpp
g++ -c `pkg-config --cflags --libs ginac` lev4.cpp
g++ -c `pkg-config --cflags --libs ginac` lev5.cpp
g++ `pkg-config --cflags --libs ginac` -o test.exe test.cpp trm2F.o my_fns.o trmchk.o lev1.o lev2.o lev3.o lev4.o lev5.o

