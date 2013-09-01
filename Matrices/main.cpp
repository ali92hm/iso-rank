

#include "DenseMatrix1D.h"
#include "DenseMatrix2D.h"
#include "SymMatrix.h"
#include <vector>
#include <sstream>
#include <string>
#include <ctime>
#ifdef USE_MPI
#include "mpi.h"
#endif

std::string G_DIR_PATH = "../Sample input/";
std::string G_FILE_EXTENSION = ".dat";
int G_NUMBER_OF_FILES = 1;

#if INT_MATRIX
typedef int DataType;
#elif DOUBLE_MATRIX
typedef double DataType;
#elif LONG_DOUBLE
typedef long double DataType;
#else
typedef float DataType;
#endif

template <typename T>
bool equal(SymMatrix<T> mat1, DenseMatrix2D<T> mat2)
{
	if (mat1.getNumberOfRows() != mat2.getNumberOfRows() || mat1.getNumberOfColumns() != mat2.getNumberOfColumns())
	{
		return false;
	}

	for (int i = 0; i < mat1.getNumberOfRows(); i++)
	{
		for (int j = 0; j < mat1.getNumberOfColumns(); j++)
		{
			if (mat1(i,j) != mat2(i,j))
			{
				return false;
			}
		}
	}
	return true;
}

template <typename T>
bool equal(SymMatrix<T> mat1, DenseMatrix1D<T> mat2)
{
	if (mat1.getNumberOfRows() != mat2.getNumberOfRows() || mat1.getNumberOfColumns() != mat2.getNumberOfColumns())
	{
		return false;
	}

	for (int i = 0; i < mat1.getNumberOfRows(); i++)
	{
		for (int j = 0; j < mat1.getNumberOfColumns(); j++)
		{
			if (mat1(i,j) != mat2(i,j))
			{
				return false;
			}
		}
	}
	return true;
}

template <typename T>
bool equal(DenseMatrix2D<T> mat1, DenseMatrix1D<T> mat2)
{
	if (mat1.getNumberOfRows() != mat2.getNumberOfRows() || mat1.getNumberOfColumns() != mat2.getNumberOfColumns())
	{
		return false;
	}

	for (int i = 0; i < mat1.getNumberOfRows(); i++)
	{
		for (int j = 0; j < mat1.getNumberOfColumns(); j++)
		{
			if (mat1(i,j) != mat2(i,j))
			{
				return false;
			}
		}
	}
	return true;
}

template<typename T>
bool qeual(std::vector<T> vec1, std::vector<T> vec2)
{
	if (vec1.size() != vec2.size())
	{
		return false;
	}

	for(int i = 0; i < vec1.size(); i++)
	{
		if (vec1[i] != vec2[i])
		{
			return false;
		}
	}
	return true;
}

template<typename T>
bool qeual(std::vector<SparseElement<T> > vec1, std::vector<SparseElement<T> > vec2)
{
	if (vec1.size() != vec2.size())
	{
		return false;
	}

	for(int i = 0; i < vec1.size(); i++)
	{
		if (vec1[i] == vec2[i])
		{
			
		}
		else
		{
			return false;
		}
	}
	return true;
}

#ifdef USE_MPI
int main(int argc, char *argv[])
{
	 int num_procs;
    int ID;
 	MPI_Status stat;
    
 	/*
 	 * MPI constant Tags
 	 */
 	const int MASTER_ID = 0;
    const int TAG_1 = 4;
    const int TAG_2 = 10;
    const int TAG_3 = 15;
    /*
     * MPI Initialization calls
     */
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        std::cout << "Failed To Initialize MPI" << std::endl;
        //MPI_Abort();
    }
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank (MPI_COMM_WORLD, &ID);
	std::vector<DenseMatrix1D<DataType>* >input_graphs;



		std::ostringstream itos_converter;

		/*
		 * Reading the graphs and storing them
		 */
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
				input_graphs.push_back(new DenseMatrix1D<DataType>(itos_converter.str()));
				itos_converter.str(""); //clearing the stream
				itos_converter.clear();
			}
			catch (std::exception& e)
			{
				std::cerr <<"Exception: " << e.what() << '\n' << std::endl;
				itos_converter.str("");
				itos_converter.clear();
			}
		}

	if (ID == 0)
	{
		for(int i = 0; i < input_graphs.size(); i++)
		{
			input_graphs[i]->MPI_Send_Matrix(1,4);
			// input_graphs[i]->MPI_Bcast_Send_Matrix(0);

		}

		for(int i = 0; i < input_graphs.size(); i++)
		{
			std::cout <<  "Master\n" << *input_graphs[i] << std::endl;
		}
	}
	else
	{
		std::vector<DenseMatrix2D<DataType>* >recv_graphs;
		for(int i = 0; i < 1; i++)
		{
			recv_graphs.push_back(new DenseMatrix2D<DataType>(0,4,stat));
			// recv_graphs.push_back(new DenseMatrix2D<DataType>(0,stat));
		}

		for(int i = 0; i < recv_graphs.size(); i++)
		{
			std::cout <<  "Slave\n" << *recv_graphs[i] << std::endl;
		}

		for(int i = 0; i < recv_graphs.size(); i++)
		{
			std::cout <<  equal(*recv_graphs[i] ,*input_graphs[i])  << std::endl;
		}

	}



	MPI_Finalize();
	return 0;
}
#else
int main(int argc, char *argv[])
{
	    std::vector<SymMatrix<DataType> >input_graphs_1D;
	    std::vector<SymMatrix<DataType> >input_graphs_2D;
	    std::ostringstream itos_converter;
    	/*
    	 * Reading the graphs and storing them
    	 */
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
				input_graphs_1D.push_back(SymMatrix<DataType>(itos_converter.str()));
				input_graphs_2D.push_back(SymMatrix<DataType>(itos_converter.str()));
				itos_converter.str(""); //clearing the stream
				itos_converter.clear();
			}
			catch (std::exception& e)
			{
				std::cerr <<"Exception: " << e.what() << '\n' << std::endl;
				itos_converter.str("");
				itos_converter.clear();
			}
		}


		
		std::cout << "Raw matrices are equal: \n"<< qeual(input_graphs_1D[0], input_graphs_2D[0]) << std::endl;
		// std::cout << input_graphs_1D[0].isSquare() << std::endl;
		// std::cout << input_graphs_1D[0].isSymmetric()<< std::endl;
		// std::cout << input_graphs_1D[0].getNumberOfRows()<< std::endl;
		// std::cout << input_graphs_1D[0].getNumberOfColumns()<< std::endl;
		std::cout << "Sparse form are equal: \n" << qeual(input_graphs_1D[0].getSparseForm(), input_graphs_2D[0].getSparseForm()) << std::endl;
		// std::cout << "FrobNorm are qeual: \n" << (input_graphs_1D[0].getFrobNorm() == input_graphs_2D[0].getFrobNorm()) << std::endl;
		// std::cout << "Transpose are equal: \n"<< qeual(input_graphs_1D[0].transpose(), input_graphs_2D[0].transpose()) << std::endl;
		// std::cout << "Kron Products are equal: \n"<< qeual(input_graphs_1D[0].kron(input_graphs_1D[0]), input_graphs_2D[0].kron(input_graphs_2D[0])) << std::endl;
		// std::cout << "Scattered Selections are qeual \n" << qeual(input_graphs_1D[0].getScatteredSelection(vec_A,vec_A), input_graphs_2D[0].getScatteredSelection(vec_A,vec_A)) << std::endl;
		// std::cout << "Sum of rows are equal \n" << qeual(input_graphs_1D[0].getSumOfRows(), input_graphs_2D[0].getSumOfRows()) << std::endl;
		// std::cout << "Neighbors \n" << qeual(input_graphs_1D[0].getNeighbors(4), input_graphs_2D[0].getNeighbors(4)) << std::endl;
		// std::cout << "diagonalVectorTimesMatrix \n" << qeual(input_graphs_1D[0].diagonalVectorTimesMatrix(vec_B), input_graphs_2D[0].diagonalVectorTimesMatrix(vec_B)) << std::endl;
		// std::cout << "matrixTimesDiagonalVector \n" << qeual(input_graphs_1D[0].matrixTimesDiagonalVector(vec_B), input_graphs_2D[0].matrixTimesDiagonalVector(vec_B)) << std::endl;
		// std::cout << "+ \n" << qeual(input_graphs_1D[0] + input_graphs_1D[0], input_graphs_2D[0] + input_graphs_2D[0]) << std::endl;
		// std::cout << "- \n" << qeual(input_graphs_1D[0] - input_graphs_1D[0], input_graphs_2D[0] - input_graphs_2D[0]) << std::endl;
		// std::cout << "* \n" << qeual(input_graphs_1D[0] * input_graphs_1D[0], input_graphs_2D[0] * input_graphs_2D[0]) << std::endl;

		// for(int i = 0; i < input_graphs_2D.size(); i++)
		// {
		// 	std::cout << input_graphs_1D[i] << std::endl;
		// 	std::cout << input_graphs_2D[i] << std::endl;

		// }
	return 0;
}
#endif