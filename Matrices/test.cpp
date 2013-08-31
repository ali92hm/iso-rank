

#include "DenseMatrix1D.h"
// #include "DenseMatrix1D.h"
// #include "DenseMatrix2D.h"
#include <vector>
#include <sstream>
#include <ostream>
#include <string>
#include <ctime>
#include "mpi.h"

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
std::ostream& operator<<(std::ostream& stream, std::vector<T>& vec)
{
	stream << "[";
	for (int i = 0; i < vec.size() - 1; i++)
	{
		stream << vec[i] << ", ";
	}
	stream << vec[vec.size() - 1] << "]" ;
	return stream;
}

template <typename T>
void printArray(T* array, int size)
{
	std::cout << "[";
	for (int i = 0; i < size - 1; i++)
	{
		std::cout << array[i] << ", ";
	}
	std::cout << array[size - 1] << "]" << std::endl;
}


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
    // std::cout << sizeof(SparseElement<float>) << std::endl;
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
	

	if (ID == 0)
	{
		std::vector<SparseElement<DataType> > sparse_form(5);
		for (int i=0; i < sparse_form.size(); i++)
		{
			sparse_form[i] = SparseElement<DataType>(rand()%10, rand()%10, rand()%10);
		}
		std::cout << "Master: " << sparse_form << std::endl;
		MPI_Send(&sparse_form[0], sparse_form.size() * sizeof(SparseElement<DataType>), MPI_BYTE, 1, 5, MPI_COMM_WORLD);
		std::cout << "Object was sent" << std::endl;

	}
	else
	{
		std::vector<SparseElement<DataType> > sparse_form(5);
		std::cout << "slave: " << sparse_form << std::endl;
		MPI_Recv(&sparse_form[0], 5 * sizeof(SparseElement<DataType>), MPI_BYTE, 0, 5, MPI_COMM_WORLD, &stat);
		std::cout << "slave: " << sparse_form << std::endl;
	}
	while(true)
	{
		
	}

	MPI_Finalize();
	return 0;
}
