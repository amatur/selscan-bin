CC=g++
# CFLAGS=-w -std=c++11 -O3  -ftree-vectorize -pthread
CFLAGS= -O3 -pthread
#-g
SRC_DIR := ./src
OBJ_DIR := ./
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
LDFLAGS := -pthread

#g++ src/main.cpp -O3 -std=c++0x -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=1 -o main

#-w suppresses warning
#enable -DDEBUGMODE for debugging

all: make_directories selbin

selbin: $(OBJ_FILES) 
	g++ $(LDFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp src/CLI11.hpp
	g++ $(CFLAGS) -c -o $@ $<
# selbin: 
# 	$(CC) $(CFLAGS) -o selbin  src/main.cpp src/nsl.cpp


# $(CC) $(CFLAGS) -c src/*.cpp
# $(CC) *.o -o selbin
# $(CC) $(CFLAGS) -o selbin  src/main.cpp src/nsl.cpp

.PHONY: make_directories
make_directories:
	mkdir -p ../bin/

clean:
	rm -rf *.o selbin