//
//  main.cpp
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/11/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#include "SparseMatrix.h"
#include "IsoRank.h"
#include <ctime>

/*
 * generate a random square matrix of type int
 * @pram: int size of the matrix
 */
template <typename DT>
void getRandomMatrix(SparseMatrix<DT>& matrix )
{
    std::srand(1);
    for(int i=0; i< matrix.getNumberOfRows() ; i++ )
    {
        for(int j=0; j < i; j++)
        {
            matrix[i][j] = matrix[j][i] = std::rand() % 2;
        }
    }
}



int main(int argc, const char * argv[])
{
//    std::string DIR_PATH = "/Users/AliHM/Documents/Course Material/Summer 13 REU/graphs/";
//    std::string EXTENSION = ".dat";
//    int NUMBER_OF_FILES = 1;
//    std::string File = "1";

    SparseMatrix<float> matrix_A(3,3);
    getRandomMatrix(matrix_A);
    std::cout<< matrix_A << std::endl;
    SparseMatrix<float> matrix_B(5,5);
    getRandomMatrix(matrix_B);
    std::cout<< matrix_B << std::endl;
    isoRank(matrix_A, matrix_B);

    
    
    return 0;
}




