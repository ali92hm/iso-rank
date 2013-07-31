CC=mpicxx
CC_SEQ=g++
CFLAGS= -g -m32 -DUSE_MPI  -DARPACK  
CFLAGS_SEQ= -g -m32  -DARPACK 
ARPACK_DIR= $(HOME)/reu_share/lib/arpack++/
INCLUDE= -I$(ARPACK_DIR)/include/ -I$(ARPACK_DIR)/examples/matrices/nonsym -I$(ARPACK_DIR)/examples/matrices/sym -I/usr/local/include/eigen3/
LIBRARIES= /share/apps/lib/libarpack.a /share/apps/lib/libsuperlu_4.3.a /usr/lib/libblas.so.3.2.1 /usr/lib/liblapack.so.3.2.1 /share/apps/lib/libf2c.a -lm

all:
	$(CC) $(CFLAGS) $(INCLUDE) $(INCLUDE_MPI) -c main.cpp 
	$(CC) $(CFLAGS) $(INCLUDE) -o IsoRank main.o $(LIB_MPI) $(LIBRARIES)
	$(CC_SEQ) $(CFLAGS_SEQ) $(INCLUDE) -c main_seq.cpp 
	$(CC_SEQ) $(CFLAGS_SEQ) $(INCLUDE) -o IsoRank_seq main_seq.o $(LIBRARIES)
	