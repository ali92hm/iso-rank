

// #include "SymMatrix.h"
#include "DenseMatrix1D.h"
// #include "DenseMatrix2D.h"
#include <vector>
#include <sstream>
#include <string>
#include <ctime>
#ifdef USE_MPI
#include "mpi.h"
#endif
#ifdef __linux__
std::string G_DIR_PATH = "/home/ali/lib/input/";
#elif defined __APPLE__
std::string G_DIR_PATH = "/Users/AliHM/Documents/Course Material/Summer 13 REU/graphs/";
#endif

std::string G_FILE_EXTENSION = ".dat";
int G_NUMBER_OF_FILES = 2;

#if INT_MATRIX
typedef int DataType;
#elif DOUBLE_MATRIX
typedef double DataType;
#elif LONG_DOUBLE
typedef long double DataType;
#else
typedef float DataType;
#endif

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

	for(int i = 0; i < input_graphs.size(); i++)
	{
		std::cout << *input_graphs[i] << std::endl;
	}
	MPI_Finalize();
	return 0;
}
#else
int main(int argc, char *argv[])
{
	    std::vector<SymMatrix<DataType>* >input_graphs;
	    std::ostringstream itos_converter;
    	/*
    	 * Reading the graphs and storing them
    	 */
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
				input_graphs.push_back(new SymMatrix<DataType>(itos_converter.str()));
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
		for(int i = 0; i < input_graphs.size(); i++)
		{
			std::cout << *input_graphs[i] << std::endl;
		}
	return 0;
}
#endif