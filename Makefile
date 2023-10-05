LD 				 := g++
CXX				 := c++
CXX_FLAGS 		 := -std=c++17 -ggdb 
CXX_STATIC_FLAGS := -std=c++17 -ggdb -static
CXX_PREPROCESSOR := -DCOMPILE_AS_PROGRAM

BIN     		:= bin
SRC     		:= source
INCLUDE 		:= include
PY_INCLUDES		:= $(shell python3-config --includes)
INTERFACEDIR 	:= python


LIBRARIES   		:= -lgsl -lgslcblas -fopenmp -lfftw3
EXECUTABLE  		:= main
STATICLIBRARIES 	:= -ldl
STATICEXECUTABLE 	:= smain
SWIGINTERFACE		:= icegiant.i


all: $(BIN)/$(EXECUTABLE)

run: clean all
		clear
		./$(BIN)/$(EXECUTABLE)

static:
	$(CXX) $(CXX_STATIC_FLAGS) $(CXX_PREPROCESSOR) $(SRC)/*.cpp -I$(INCLUDE) -o $(BIN)/$(STATICEXECUTABLE) $(LIBRARIES) $(STATICLIBRARIES)


$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp 
		$(CXX) $(CXX_FLAGS) $(CXX_PREPROCESSOR) -I$(INCLUDE) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*

.python:
	swig -python -c++ $(INTERFACEDIR)/$(SWIGINTERFACE)
	$(CXX) -c $(CXX_FLAGS) $(PY_INCLUDES) -I$(INCLUDE) $(SRC)/Binary.cpp $(SRC)/LISA.cpp $(SRC)/utils.cpp $(INTERFACEDIR)/icegiant_wrap.cxx $(LIBRARIES) -fPIC
	$(LD) -shared *.o -o _icegiant.so
	-rm *.o
	-mv _icegiant.so $(INTERFACEDIR)

