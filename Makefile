all:vertex IsoRank

vertex: vertex.h vertex.cpp
	g++ -Wall -pedantic -g vertex.cpp -o vertex

IsoRank: IsoRank.h vertex.h SparseMatrix.h Tarjan.h
	g++ -Wall -pedantic -g  vertex main.cpp -o IsoRank
