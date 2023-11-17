CC=g++
# CFLAGS=-w -std=c++11 -O3  -ftree-vectorize -pthread
CFLAGS= -Ofast -pthread -m64 -mmmx -msse -msse2 -ftree-vectorize -fopenmp
#-g
SRC_DIR := ./src
OBJ_DIR := ./
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
LDFLAGS := -pthread


#for osx systems
CC = g++
G++FLAG = -Ofast -m64 -mmmx -msse -msse2 -fopenmp
LINK_OPTS2 = -lpthread -lz -fopenmp
I_PATH = -I../include
L_PATH = ../lib/osx

#for linux systems
#CC = g++
#G++FLAG = -O3 -m64 -mmmx -msse -msse2
#LINK_OPTS2 = -pthread -lz
#I_PATH = -I../include
#L_PATH = ../lib/linux

#for windows, using MinGW build environment
#CC = g++.exe
#G++FLAG = -DPTW32_STATIC_LIB -O3 -static-libgcc -static-libstdc++
#LINK_OPTS2 = ../lib/win32/libpthreadGC2.a ../lib/win32/libz.a 
#I_PATH = -I../include -I../include/win32 
#L_PATH = ../lib/win32

#For static linking of norm program to gsl libs
LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a
#g++ src/main.cpp -O3 -std=c++0x -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=1 -o main

#-w suppresses warning
#enable -DDEBUGMODE for debugging

all: make_directories selbin

selbin: $(OBJ_FILES) 
	g++ $(LDFLAGS) -o $@ $^ $(LINK_OPTS2)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp src/CLI11.hpp
	g++ $(CFLAGS) -c -o $@ $< $(I_PATH)
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
