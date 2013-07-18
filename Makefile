CC=mpicxx
CFLAGS= -g -DEIGEN
INCLUDE= -I/usr/local/include/eigen3/

all:
	$(CC) $(CFLAGS) $(INCLUDE) $(INCLUDE_MPI) -c main.cpp 
	$(CC) $(CFLAGS) $(INCLUDE) -o main main.o $(LIB_MPI) $(LIBRARIES)
