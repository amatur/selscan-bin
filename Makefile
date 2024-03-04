CC=g++
CFLAGS= -Ofast -pthread -m64  -ftree-vectorize -fopenmp -mmmx -msse -msse2
SRC_DIR := ./src
OBJ_DIR := ./
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
LDFLAGS := -pthread

LINK_OPTS2 = -lpthread -lz -fopenmp
I_PATH = -I../include
L_PATH = ../lib/linux

all: make_directories selbin

selbin: $(OBJ_FILES) 
	$(CC) $(LDFLAGS) -o $@ $^ $(LINK_OPTS2)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp src/CLI11.hpp
	$(CC) $(CFLAGS) -c -o $@ $< $(I_PATH)
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
