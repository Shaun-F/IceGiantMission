CXX	:= c++
CXX_FLAGS := -std=c++17 -ggdb
CXX_STATIC_FLAGS := -std=c++17 -ggdb -static

BIN     := bin
SRC     := source
INCLUDE := include

LIBRARIES   		:= -lgsl -lgslcblas -fopenmp -lfftw3
EXECUTABLE  		:= main
STATICLIBRARIES 	:= -ldl
STATICEXECUTABLE 	:= smain


all: $(BIN)/$(EXECUTABLE)

run: clean all
		clear
		./$(BIN)/$(EXECUTABLE)

static:
	$(CXX) $(CXX_STATIC_FLAGS) $(SRC)/*.cpp -I$(INCLUDE) -o $(BIN)/$(STATICEXECUTABLE) $(LIBRARIES) $(STATICLIBRARIES)


$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp 
		$(CXX) $(CXX_FLAGS) -I$(INCLUDE) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*
