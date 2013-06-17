all:IsoRank


IsoRank: IsoRank.h vertex.h SparseMatrix.h Tarjan.h
	g++ -Wall -pedantic -g  main.cpp -o IsoRank
