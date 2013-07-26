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
int G_NUMBER_OF_FILES = 5000;

typedef int DataType;
double timeElapsed(std::clock_t start, std::clock_t end)
{
	return (double) (end - start) / CLOCKS_PER_SEC * 1000.0;
}

using namespace std;
int main(int argc, char *argv[])
{
	int rank = atoi(argv[1]);
	int num_procs;
	std::clock_t time_start;
    std::clock_t time_end;

	// if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
 //    {
 //        std::cout << "Failed To Initialize MPI" << std::endl;
 //        //MPI_Abort();
 //    }
 //    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
 //    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		std::ostringstream itos_converter;
		std::vector<SymSparseMatrix<float> > input_graphs;
		time_start = std::clock();
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
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
		time_end = std::clock();
		std::cout << "Rank " << rank <<": " <<  input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
			<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;

		time_start = std::clock();
		for (int i =0 ; i < input_graphs.size(); i++)
		{
			for (int j = i ; j < input_graphs.size(); j++)
			{
				SymSparseMatrix<float> new_rep = input_graphs[i].kron(input_graphs[j]);
			}
		}

		time_end = std::clock();
		std::cout << "Rank " << rank <<": Computed Kron prod successfully for " << G_NUMBER_OF_FILES << " graphs in "
												<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;
	}

	if (rank == 1)
	{
		
		std::ostringstream itos_converter;
		std::vector<SparseMatrix<float> > input_graphs;
		time_start = std::clock();
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
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
		time_end = std::clock();
		std::cout << "Rank " << rank <<": " <<  input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
			<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;

		time_start = std::clock();	
		for (int i = 0 ; i < input_graphs.size(); i++)
		{
			for (int j = i ; j < input_graphs.size(); j++)
			{
			  	SparseMatrix<float> new_rep = input_graphs[i].kron(input_graphs[j]);
			}
		}

		time_end = std::clock();
		std::cout << "Rank " << rank <<": Computed Kron prod successfully for " << G_NUMBER_OF_FILES << " graphs in "
												<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;
	}

	if (rank == 2)
	{
		
		std::ostringstream itos_converter;
		std::vector<Matrix<float> > input_graphs;
		time_start = std::clock();
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
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
		time_end = std::clock();
		std::cout << "Rank " << rank <<": " <<  input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
			<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;

		time_start = std::clock();
		for (int i = 0 ; i < input_graphs.size(); i++)
		{
			for (int j = i ; j < input_graphs.size(); j++)
			{
			  	Matrix<float> new_rep = input_graphs[i].kron(input_graphs[j]);
			}
		}
		time_end = std::clock();
		std::cout << "Rank " << rank <<": Computed Kron prod successfully for " << G_NUMBER_OF_FILES << " graphs in "
												<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;
	}

	if (rank == 3)
	{
		
		std::ostringstream itos_converter;
		std::vector<SymMatrix<float> > input_graphs;
		time_start = std::clock();
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
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
		time_end = std::clock();
		std::cout << "Rank " << rank <<": " <<  input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
			<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;

		time_start = std::clock();
		for (int i = 0 ; i < input_graphs.size(); i++)
		{
			for (int j = i ; j < input_graphs.size(); j++)
			{
			  	SymMatrix<float> new_rep = input_graphs[i].kron(input_graphs[j]);

			}
		}
		time_end = std::clock();
		std::cout << "Rank " << rank <<": " <<  input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
			<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;
	}

	if (rank == 4)
	{
		std::ostringstream itos_converter;
		std::vector<DenseMatrix<float> > input_graphs_2D;
		time_start = std::clock();
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
				input_graphs_2D.push_back(DenseMatrix<float>(itos_converter.str()));
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
		time_end = std::clock();
		std::cout << "Rank " << rank <<": " <<  input_graphs_2D.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
			<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;

		time_start = std::clock();
		for (int i = 0 ; i < input_graphs_2D.size(); i++)
		{
			for (int j = i ; j < input_graphs_2D.size(); j++)
			{
				DenseMatrix<float> old_rep = input_graphs_2D[i].kron(input_graphs_2D[j]);
			}
		}
		time_end = std::clock();
		std::cout << "Rank " << rank <<": " <<  input_graphs_2D.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
			<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;
	}

	std::cout << "Process: " << rank << " terminated." << endl;
	// MPI_Finalize();
	return 0;
}
