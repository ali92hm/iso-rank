CC=g++
CFLAGS= -m32  -g
ARPACK_DIR= $(HOME)/reu_share/lib/arpack++
MPI_LIB_DIR=/usr/lib/openmpi/lib
INCLUDE= -I$(ARPACK_DIR)/include/ -I$(ARPACK_DIR)/examples/matrices/nonsym -I$(ARPACK_DIR)/examples/matrices/sym -I/usr/include/openmpi-x86_64/
INCLUDE_MPI= -pthread -L$(MPI_LIB_DIR) -lmpi_cxx -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl
LIBRARIES= /share/apps/lib/libarpack.a /share/apps/lib/libsuperlu_4.3.a /usr/lib/libblas.so.3.2.1 /usr/lib/liblapack.so.3.2.1 /share/apps/lib/libf2c.a -lm 
LIB_MPI= $(MPI_LIB_DIR)/libmpi_cxx.so.1 $(MPI_LIB_DIR)/libmpi_cxx.so.1.0.1 $(MPI_LIB_DIR)/libmpi.so.1 $(MPI_LIB_DIR)/libmpi.so.1.0.1 $(MPI_LIB_DIR)/libopen-rte.so.2 $(MPI_LIB_DIR)/libopen-rte.so.2.0.0 $(MPI_LIB_DIR)/libopen-pal.so.2.0.0 $(MPI_LIB_DIR)/libopen-pal.so.2.0.0

all:
	$(CC) $(CFLAGS) $(INCLUDE) $(INCLUDE_MPI) -c main.cpp 
	$(CC) $(CFLAGS) $(INCLUDE) -o main main.o $(LIB_MPI) $(LIBRARIES)
