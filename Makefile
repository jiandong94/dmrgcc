## $@ target
## $^ all prerequisites
## $< first prerequisite

#path
ROOT_DIR = $(shell pwd)
BIN_DIR = $(ROOT_DIR)/bin
SRC_DIR = $(ROOT_DIR)/src
INC_DIR = $(ROOT_DIR)/include
OBJ_DIR = $(ROOT_DIR)/obj

SUBDIRS=$(shell ls -l | grep ^d | awk '{if($$9 == "src") print $$9}')

CC = g++
CFLAGS = -fopenmp -m64 -std=c++11 -fPIC -O2 -msse2 -msse3 -msse4 
#CFLAGS = -fopenmp -std=c++11 -m64 -fPIC -O3
#CFLAGS = -fopenmp -m64 -std=c++11 -fPIC -g
MKL_INCLUDE_PATH = /opt/intel/mkl/include 
MKL_LIBRARY_PATH = /opt/intel/mkl/lib/intel64 -Wl,--start-group -lmkl_gnu_thread -lmkl_core -lmkl_intel_lp64 -Wl,--end-group
INCLUDE = -I$(INC_DIR) -I$(SRC_DIR) -I$(MKL_INCLUDE_PATH)
LIBRARY = -L$(MKL_LIBRARY_PATH)
vpath *.o $(OBJ_DIR)

export BIN_DIR SRC_DIR INC_DIR OBJ_DIR CC CFLAGS INCLUDE LIBRARY


all:$(SUBDIRS) TEST MODEL

$(SUBDIRS):ECHO
	make -C $@
TEST:ECHO
	make -C test
MODEL:ECHO
	make -C model
ECHO:
	@echo $(SUBDIRS)

clean:
	@rm -rf $(OBJ_DIR)/*.o
	@rm -rf $(BIN_DIR)/*
