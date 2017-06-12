
CXX = c++
CXXFLAGS = -pthread -std=c++14
OBJS = src/args.o src/dictionary.o src/productquantizer.o src/matrix.o src/qmatrix.o src/vector.o src/model.o src/utils.o src/fasttext.o 
INCLUDES = -I./src

opt: CXXFLAGS += -O3 -funroll-loops
opt: ARGfast

debug: CXXFLAGS += -g -O0 -fno-inline
debug: ARGfast

args.o: src/args.cc src/args.h
	$(CXX) $(CXXFLAGS) -c src/args.cc -o src/args.o

dictionary.o: src/dictionary.cc src/dictionary.h src/args.h
	$(CXX) $(CXXFLAGS) -c src/dictionary.cc -o src/dictionary.o

productquantizer.o: src/productquantizer.cc src/productquantizer.h src/utils.h
	$(CXX) $(CXXFLAGS) -c src/productquantizer.cc -o src/productquantizer.o

matrix.o: src/matrix.cc src/matrix.h src/utils.h
	$(CXX) $(CXXFLAGS) -c src/matrix.cc -o src/matrix.o

qmatrix.o: src/qmatrix.cc src/qmatrix.h src/utils.h
	$(CXX) $(CXXFLAGS) -c src/qmatrix.cc -o src/qmatrix.o

vector.o: src/vector.cc src/vector.h src/utils.h
	$(CXX) $(CXXFLAGS) -c src/vector.cc -o src/vector.o

model.o: src/model.cc src/model.h src/args.h
	$(CXX) $(CXXFLAGS) -c src/model.cc -o src/model.o

utils.o: src/utils.cc src/utils.h
	$(CXX) $(CXXFLAGS) -c src/utils.cc -o src/utils.o

fasttext.o: src/fasttext.cc src/*.h
	$(CXX) $(CXXFLAGS) -c src/fasttext.cc -o src/fasttext.o

ARGfast: $(OBJS) src/map.cpp
	$(CXX) $(CXXFLAGS) $(OBJS) $(INCLUDES) src/map.cpp -o bin/predX 

clean:
	rm -rf src/*.o bin/predX
