CC=mpicxx
CFLAGS= -Wall -g 
LDFLAGS = -lm


point.o : Point.cpp Point.h
	${CC} ${CFLAGS} -c Point.cpp

main.o : test.c Point.h
	${CC} ${CFLAGS} -c test.cpp 

IsoRank_seq: main.o point.o
	${CC} ${CFLAGS} point.o main.o ${LDFLAGS} -o IsoRank_seq
	