CC=g++
CFLAGS= -m32 -ggdb
ARPACK_DIR= $(HOME)/arpack++/arpack++
INCLUDE= -I$(ARPACK_DIR)/include/ -I$(ARPACK_DIR)/examples/matrices/nonsym -I$(ARPACK_DIR)/examples/matrices/sym
LIBRARIES= /share/apps/lib/libarpack.a /share/apps/lib/libsuperlu_4.3.a /usr/lib/libblas.so.3.2.1 /usr/lib/liblapack.so.3.2.1 /share/apps/lib/libf2c.a -lm


all:
	$(CC) $(CFLAGS) $(INCLUDE) -c main.cpp
	$(CC) $(CFLAGS) -o main main.o $(LIBRARIES)
