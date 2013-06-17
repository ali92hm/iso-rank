cmake_minimum_required(VERSION 2.6)

project (IsoRank) 

SET(sampleName IsoRank)

SET(CMAKE_C_FLAGS "-m32 -ggdb -D_DEBUG -Wextra -Wall -Wunreachable-code -O0 -D_LARGEFILE64_SOURCE")

add_executable (IsoRank main.cpp IsoRank.h util.h vertex.cpp vertex.h SparseMatrix.h Tarjan.h)
