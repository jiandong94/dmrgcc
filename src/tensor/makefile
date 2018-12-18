CCCOM = g++ -fopenmp -m64 -std=c++11 -fPIC
CFLAGS = -O4 -msse2 -msse3 -msse4 
INCLUDES = -I/home/chenjd/Code/GIT/dmrgcc/src/ 
MKL_INCLUDES = -I//opt/intel/mkl/include 
MKL_LIBRARY = -L/opt/intel/mkl/lib/intel64 -Wl,--start-group -lmkl_gnu_thread -lmkl_core -lmkl_intel_lp64 -Wl,--end-group
#MKL_LIBRARY = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_rt -lmkl_core -liomp5 -lpthread
#-lmkl_sequential_thread
#-I/home/chenjd/Code/GIT/DMRGByWu/BETA/src 


#VPATH = /home/chenjd/Code/GIT/DMRGByWu/BETA/src/Tensor
OUTNAME = test_real_matrix_block
OBJECTS = test_real_matrix_block.o real_matrix.o real_matrix_block.o
#OUTNAME = test_real_matrix
#OBJECTS = test_real_matrix.o real_matrix.o real_matrix_block.o

all : ${OBJECTS}
	${CCCOM} ${CFLAGS} -o ${OUTNAME} $^ ${INCLUDES} ${MKL_LIBRARY}

%.o : %.cc
	${CCCOM} ${CFLAGS} ${INCLUDES} ${MKL_INCLUDES} -c $< -o $@

clean :
	rm -rf *.o
