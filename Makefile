CC=mpicxx
# Add -DSEQ for sequential code
# Add -DNODE_PAIR for node pair method
# default method is broadcast
CFLAGS= -O3 -m32 -DARPACK -DUSE_MPI
ARPACK_DIR= $(HOME)/reu_share/lib/arpack++/
INCLUDE= -I$(ARPACK_DIR)/include/ -I$(ARPACK_DIR)/examples/matrices/nonsym -I$(ARPACK_DIR)/examples/matrices/sym -I/usr/local/include/eigen3/
LIBRARIES= /share/apps/lib/libarpack.a /share/apps/lib/libsuperlu_4.3.a /usr/lib/libblas.so.3.2.1 /usr/lib/liblapack.so.3.2.1 /share/apps/lib/libf2c.a -lm

all:
	$(CC) $(CFLAGS) $(INCLUDE) $(INCLUDE_MPI) -c Vertex.cpp 
	$(CC) $(CFLAGS) $(INCLUDE) $(INCLUDE_MPI) -c main.cpp 
	$(CC) $(CFLAGS) $(INCLUDE) -o IsoRank main.o Vertex.o $(LIB_MPI) $(LIBRARIES)
	rm -f *.o

	