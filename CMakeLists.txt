cmake_minimum_required(VERSION 2.6)

project (IsoRank) 

SET(sampleName IsoRank)

add_executable (IsoRank main.cpp IsoRank.h util.h vertex.cpp vertex.h SparseMatrix.h Tarjan.h)
