#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "SymSparseMatrix.h"
#include "Matrix2D.h"
#include <time.h>
#include "../Tarjan.h"

// #include "SparseMatrix.h"
// #include "SymSparseMatrix.h"

#ifdef __linux__
std::string G_DIR_PATH = "/home/ali/workspace/ex/input/";
#elif defined __APPLE__
std::string G_DIR_PATH = "/Users/AliHM/Documents/Course Material/Summer 13 REU/graphs/";
#endif
std::string G_FILE_EXTENSION = ".dat";
int G_NUMBER_OF_FILES = 50;

typedef int DataType;

template <typename DT>
bool equal(SparseMatrix<DT>& mat1, SymSparseMatrix<DT>& mat2)
{
	for(int i=0; i < mat2.getSize(); i++)
	{
		for(int j = 0; j < mat2.getSize(); j++)
		{
			if (mat2(i,j) != mat1[i][j])
			{
				return false;
			}
		}
	}
	return true;
}

using namespace std;
int main(int argc, char *argv[])
{
	srand (2);
	std::ostringstream itos_converter;
	std::vector<SparseMatrix<int> > input_graphs_2D;
	std::vector<SymSparseMatrix<int> > input_graphs;


	for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
	{
		try
		{
			itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
			input_graphs_2D.push_back(SparseMatrix<int>(itos_converter.str()));
			input_graphs.push_back(SymSparseMatrix<int>(itos_converter.str()));
			//clearing the stream
			itos_converter.str("");
			itos_converter.clear();
		}
		catch (std::exception& e)
		{
			std::cerr <<"Exception: " << e.what() << '\n' << std::endl;
			itos_converter.str("");
			itos_converter.clear();
		}
	}


	for (int i =0 ; i < input_graphs.size(); i++)
	{
		for (int j = i ; j < input_graphs.size(); j++)
		{
			SparseMatrix<int>* old_rep = input_graphs_2D[i].kron(input_graphs_2D[j]);
		  	SymSparseMatrix<int> new_rep = input_graphs[i].kron(input_graphs[j]);
		  	if (!equal(*old_rep, new_rep))
		  	{
				std::cout << "Kron product is DIFFERENT "<< i  <<" " << j << std::endl;
			}

		}
	}
	return 0;
}

/*
SymSparseMatrix<int> mat1 (4);
	mat1.insert(1,2, 4);
	mat1.insert(3,0, 5);


	std::cout << mat1 << std::endl;

	for (auto it = mat1._edges.begin(); it != mat1._edges.end(); ++it )
	{
		std::cout << "First: " <<it->first << " Second" << it->second << std::endl;
		std::cout << "		" << "i: " <<it->first/mat1.getSize() << " j: " << it->first%mat1.getSize() << std::endl;
	}

*/	