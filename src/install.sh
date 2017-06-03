# FastText

cver=11

c++ -pthread -std=c++$cver -O3 -funroll-loops -c args.cc
c++ -pthread -std=c++$cver -O3 -funroll-loops -c dictionary.cc
c++ -pthread -std=c++$cver -O3 -funroll-loops -c productquantizer.cc
c++ -pthread -std=c++$cver -O3 -funroll-loops -c matrix.cc
c++ -pthread -std=c++$cver -O3 -funroll-loops -c qmatrix.cc
c++ -pthread -std=c++$cver -O3 -funroll-loops -c vector.cc
c++ -pthread -std=c++$cver -O3 -funroll-loops -c model.cc
c++ -pthread -std=c++$cver -O3 -funroll-loops -c utils.cc
c++ -pthread -std=c++$cver -O3 -funroll-loops -c fasttext.cc


c++ -pthread -std=c++$cver -O3 -funroll-loops args.o dictionary.o productquantizer.o matrix.o qmatrix.o vector.o model.o utils.o fasttext.o main.cc -o ../bin/fasttext -I$PWD

c++ -pthread -std=c++$cver -O3 -funroll-loops args.o dictionary.o productquantizer.o matrix.o qmatrix.o vector.o model.o utils.o fasttext.o map.cpp -o ../bin/ARGfast -I$PWD
