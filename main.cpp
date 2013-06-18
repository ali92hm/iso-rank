//
//  main.cpp
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/11/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#include "SparseMatrix.h"
#include "IsoRank.h"
#include <vector>
#include <ctime>




#ifdef __linux__
std::string DIR_PATH = "/export/home/reu_share/input/";
#endif

#ifdef __APPLE__
std::string DIR_PATH = "/Users/AliHM/Documents/Course Material/Summer 13 REU/graphs/";
#endif

void commandLineArgs(int argc,const char* argv[]);

int main(int argc, const char * argv[])
{

    if (argc > 1)
    {
        commandLineArgs(argc, argv);
    }
    
//    std::string EXTENSION = ".dat";
//    int NUMBER_OF_FILES = 1;
//    std::vector<SparseMatrix<float>*> input_graphs (NUMBER_OF_FILES);
//    
//    for(int i=1; i <= NUMBER_OF_FILES; i++ )
//    {
//        std::string file_name = DIR_PATH + std::to_string(NUMBER_OF_FILES) + EXTENSION;
//        input_graphs.push_back(new SparseMatrix<float> (file_name));
//    }
//    
//    for (int i=0; i < NUMBER_OF_FILES; i++)
//    {
//        for(int j=i; j< NUMBER_OF_FILES; j++ )
//        {
//            isoRank(*input_graphs[i], *input_graphs[j]);
//        }
//    }
    

    SparseMatrix<float> a (4,4);
    SparseMatrix<float> b (3,3);
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    isoRank(a, b);
    

    
    
    return 0;
}

void commandLineArgs(int argc,const char* argv[])
{
    
    
}


/*
 * generate a random square matrix of type int
 * @pram: int size of the matrix
 */
template <typename DT>
void getRandomMatrix(SparseMatrix<DT>& matrix )
{
    std::srand(3);
    for(int i=0; i< matrix.getNumberOfRows() ; i++ )
    {
        for(int j=0; j < i; j++)
        {
            matrix[i][j] = matrix[j][i] = std::rand() % 2;
        }
    }
}


