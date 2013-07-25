#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "SymSparseMatrix.h"
#include "SparseMatrix.h"
#include "Matrix2D.h"
#include "Matrix.h"
#include "SymMatrix.h"
#include <time.h>
#include "../Tarjan.h"
#include "mpi.h"

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
bool equal(DenseMatrix<DT>& mat1, SymSparseMatrix<DT>& mat2)
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

template <typename DT>
bool equal(DenseMatrix<DT>& mat1, SparseMatrix<DT>& mat2)
{
	for(int i=0; i < mat1.getNumberOfRows(); i++)
	{
		for(int j = 0; j < mat1.getNumberOfColumns(); j++)
		{

			if (mat2(i,j) != mat1[i][j])
			{
				return false;
			}
		}
	}
	return true;
}

template <typename DT>
bool equal(DenseMatrix<DT>& mat1, Matrix<DT>& mat2)
{
	for(int i = 0; i < mat1.getNumberOfRows(); i++)
	{
		for(int j = 0; j < mat1.getNumberOfColumns(); j++)
		{

			if (mat2(i,j) != mat1[i][j])
			{
				return false;
			}
		}
	}
	return true;
}

template <typename DT>
bool equal(DenseMatrix<DT>& mat1, SymMatrix<DT>& mat2)
{
	for(int i = 0; i < mat1.getNumberOfRows(); i++)
	{
		for(int j = 0; j < mat1.getNumberOfColumns(); j++)
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
	int rank;
	int num_procs;

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        std::cout << "Failed To Initialize MPI" << std::endl;
        //MPI_Abort();
    }
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		std::ostringstream itos_converter;
		std::vector<DenseMatrix<float> > input_graphs_2D;
		std::vector<SymSparseMatrix<float> > input_graphs;
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
				input_graphs_2D.push_back(DenseMatrix<float>(itos_converter.str()));
				input_graphs.push_back(SymSparseMatrix<float>(itos_converter.str()));
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
				DenseMatrix<float>* old_rep = input_graphs_2D[i].kron(input_graphs_2D[j]);
			  	SymSparseMatrix<float> new_rep = input_graphs[i].kron(input_graphs[j]);
			  	if (!equal(*old_rep, new_rep))
			  	{
					std::cout << " ID = 0: Kron product is DIFFERENT "<< i+1  <<" " << j+1 << std::endl;
				}

			}
		}
	}

	if (rank == 1)
	{
		
		std::ostringstream itos_converter;
		std::vector<DenseMatrix<float> > input_graphs_2D;
		std::vector<SparseMatrix<float> > input_graphs;
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
				input_graphs_2D.push_back(DenseMatrix<float>(itos_converter.str()));
				input_graphs.push_back(SparseMatrix<float>(itos_converter.str()));
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


		for (int i = 0 ; i < input_graphs.size(); i++)
		{
			for (int j = i ; j < input_graphs.size(); j++)
			{
				DenseMatrix<float>* old_rep = input_graphs_2D[i].kron(input_graphs_2D[j]);
			  	SparseMatrix<float> new_rep = input_graphs[i].kron(input_graphs[j]);
			  	if (!equal(*old_rep, new_rep))
			  	{
					std::cout << " ID = 1: Kron product is DIFFERENT "<< i+1  <<" " << j+1 << std::endl;
				}
			}
		}
	}

	if (rank == 2)
	{
		
		std::ostringstream itos_converter;
		std::vector<DenseMatrix<float> > input_graphs_2D;
		std::vector<Matrix<float> > input_graphs;
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
				input_graphs_2D.push_back(DenseMatrix<float>(itos_converter.str()));
				input_graphs.push_back(Matrix<float>(itos_converter.str()));
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


		for (int i = 0 ; i < input_graphs.size(); i++)
		{
			for (int j = i ; j < input_graphs.size(); j++)
			{
				DenseMatrix<float>* old_rep = input_graphs_2D[i].kron(input_graphs_2D[j]);
			  	Matrix<float> new_rep = input_graphs[i].kron(input_graphs[j]);
			  	if (!equal(*old_rep, new_rep))
			  	{
					std::cout << " ID = 2: Kron product is DIFFERENT "<< i+1  <<" " << j+1 << std::endl;
				}
			}
		}
	}

	if (rank == 3)
	{
		
		std::ostringstream itos_converter;
		std::vector<DenseMatrix<float> > input_graphs_2D;
		std::vector<SymMatrix<float> > input_graphs;
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
				input_graphs_2D.push_back(DenseMatrix<float>(itos_converter.str()));
				input_graphs.push_back(SymMatrix<float>(itos_converter.str()));
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


		for (int i = 0 ; i < input_graphs.size(); i++)
		{
			for (int j = i ; j < input_graphs.size(); j++)
			{
				DenseMatrix<float>* old_rep = input_graphs_2D[i].kron(input_graphs_2D[j]);
			  	SymMatrix<float> new_rep = input_graphs[i].kron(input_graphs[j]);
			  	if (!equal(*old_rep, new_rep))
			  	{
					std::cout << " ID = 3: Kron product is DIFFERENT "<< i+1  <<" " << j+1 << std::endl;
				}
			}
		}
	}

	std::cout << "Process: " << rank << " terminated." << endl;
	MPI_Finalize();
	return 0;
}
