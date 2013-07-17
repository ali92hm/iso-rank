//
//  main.cpp
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/11/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#include "Matricies/Matrix2D.h"
#include "IsoRank.h"
#include <vector>
#include <ctime>
#include "Tarjan.h"
#include "vertex.h"
#include "greedy_algorithms.h"
#include "mpi.h"




/*
 *
 */
#ifdef __linux__
std::string G_DIR_PATH = "/home/ali/workspace/ex/input/";
#elif defined __APPLE__
std::string G_DIR_PATH = "/Users/AliHM/Documents/Course Material/Summer 13 REU/graphs/";
#endif


int G_NUMBER_OF_FILES = 2;
std::string G_FILE_EXTENSION = ".dat";

bool G_USE_ISORANK = true;
bool G_USE_GPGM = false;
bool G_USE_GREEDY_ALG = true;
bool G_USE_CON_ENF_1 = false;
bool G_USE_CON_ENF_2 = false;
bool G_USE_CON_ENF_3 = false;
bool G_USE_CON_ENF_4 = false;


#if INT_MATRIX
typedef int DataType;
#elif DOUBLE_MATRIX
typedef double DataType;
#else
typedef float DataType;
#endif




void parseCommandLineArgs(int argc, char * argv[]);
template <typename DT>
void MPI_Send_Matrix (SparseMatrix<DT>& matrix, int dest, int tag);
SparseMatrix<DataType> MPI_Recv_Matrix (int source, int tag , MPI_Status& stat);
//void  MPI_Recv_Matrix (int source, int tag , MPI_Status& stat);



/*
 *
 */
int main(int argc, char * argv[])
{
	int num_procs;
    int rank;
 	MPI_Status stat;
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        std::cout << "fail to init" << std::endl;
        //MPI_Abort();
    }
    
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    int tag1 = 4;
    int size_tag = 3;
    /*
     *Configure the program to use the command line args
     */
    parseCommandLineArgs(argc, argv);
    std::vector<SparseMatrix<DataType> > input_graphs;
    std::clock_t time_start;
    std::clock_t time_end;
    double elapsed_time;
    /*
     * Struecutes that contain the input
     */
    if (rank == 0)
    {
    	std::cout << "Reading " << G_NUMBER_OF_FILES << " graphs from: " << G_DIR_PATH << endl;
    	time_start = std::clock();
    	ostringstream itos_converter;
		for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
		{
			try
			{
				itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
				input_graphs.push_back(SparseMatrix<DataType>(itos_converter.str()));
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
		elapsed_time = (double) (time_end - time_start) / CLOCKS_PER_SEC * 1000.0;
		std::cout << input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "<< elapsed_time << "(ms)." << endl;
		
 		int send_size = input_graphs.size();		
 		MPI_Send(&send_size, 1, MPI_INT, 1, size_tag, MPI_COMM_WORLD);
		for(int i = 0; i < input_graphs.size(); i++)
		{	
 			MPI_Send_Matrix (input_graphs[i], 1, tag1);
		}
    }
    else
    {
	  	int size_input;
     	MPI_Recv(&size_input, 1, MPI_INT, 0, size_tag, MPI_COMM_WORLD, &stat);
     	std::cout << size_input << std::endl;
    	for(int i = 0; i < size_input; i++)
		{	
			input_graphs.push_back(MPI_Recv_Matrix (0, tag1 ,stat));
			std::cout << "Recv " << input_graphs[i]  << std::endl;
		}
    	
		/*
		 * Compute n choose 2 combination of graphs
		*/
// 		time_start = std::clock();
// 		for (int i = 0; i < input_graphs.size(); i++)
// 		{
// 			for(int j = i+1; j <  input_graphs.size(); j++)
// 			{
// 				try
// 				{
// 					if (G_USE_ISORANK)
// 					{
// 						if (G_USE_GREEDY_ALG)
// 						{
// 							isoRank(input_graphs[i], input_graphs[j], 0);
// 						}
// 					
// 						else if (G_USE_CON_ENF_1)
// 						{
// 							isoRank(input_graphs[i], input_graphs[j], 1);
// 						}
// 					
// 						else if (G_USE_CON_ENF_2)
// 						{
// 							isoRank(input_graphs[i], input_graphs[j], 2);
// 						}
// 					
// 						else if (G_USE_CON_ENF_3)
// 						{
// 							isoRank(input_graphs[i], input_graphs[j], 3);
// 						}
// 						else if (G_USE_CON_ENF_4)
// 						{
// 							isoRank(input_graphs[i], input_graphs[j], 4);
// 						}
// 					}
// 				
// 					if (G_USE_GPGM)
// 					{
// 						//GPGM(input_graphs[i], input_graphs[j]);
// 					}
// 				}
// 				catch (std::exception& e)
// 				{
// 					std::cerr << "Exception: " << e.what() << std::endl;
// 				}
// 			}
// 		}
// 
// 		time_end = std::clock();
// 	    elapsed_time = (double) (time_end - time_start) / CLOCKS_PER_SEC * 1000.0;
// 		std::cout << " completed in "<< elapsed_time << "(ms)." << endl;
	}
	
	
	/*
// 	 * Deleting objects
// 	 */
// 	vector<SparseMatrix<DataType>*>::iterator it;
// 	for ( it = input_graphs.begin(); it < input_graphs.end(); ++it )
// 	{
// 		delete (*it);
// 	}
    std::cout << "Process: "<< rank << " terminated" << std::endl; 
    MPI_Finalize();
    return 0;
}

template <typename DT>
void MPI_Send_Matrix (SparseMatrix<DT>& matrix, int dest, int tag)
{
    int m = matrix.getNumberOfRows();
    int n = matrix.getNumberOfColumns();
    MPI_Send(&m, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
    MPI_Send(&n, 1, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
    MPI_Send(matrix.get1DArr(), m*n*sizeof(DT), MPI_BYTE, dest, tag + 3, MPI_COMM_WORLD);
}

SparseMatrix<DataType> MPI_Recv_Matrix (int source, int tag , MPI_Status& stat)
{
    int m, n;
    MPI_Recv(&m, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
    MPI_Recv(&n, 1, MPI_INT, source, tag + 2, MPI_COMM_WORLD, &stat);
    DataType* dataArray = new DataType[m * n];
    MPI_Recv(dataArray, m*n*sizeof(DataType), MPI_BYTE, source, tag + 3, MPI_COMM_WORLD, &stat);
    SparseMatrix<DataType> matrix (m,n, dataArray);
//     std::cout<< matrix << std::endl;
    delete [] dataArray;
    return matrix;
}

/*
 *
 */
void parseCommandLineArgs(int argc,char* argv[])
{
    std::cout << "Reading command line arguments..." << std::endl;
    
    if ( argc < 2)
    {
        std::cout <<"No arguments were passed." << std::endl;
    }
    
    for (int i = 1; i < argc; i++)
    {
        if (std::strncmp(argv[i], "-dir",4) == 0)
        {
            i++;
            G_DIR_PATH = std::string(argv[i]);
            std::cout << "Working directory was set to: " << G_DIR_PATH << std::endl;
        }
        else if (std::strncmp(argv[i], "-ext",4) == 0)
        {
            i++;
            G_FILE_EXTENSION = std::string(argv[i]);
            std::cout << "File extention was set to: " << G_FILE_EXTENSION << std::endl;
        }
        else if (std::strncmp(argv[i], "-#",2) == 0)
        {
            i++;
            int input_number = atoi(argv[i]);
            if ( input_number > 0)
            {
                G_NUMBER_OF_FILES = input_number;
                std::cout << "File extention was set to: " << G_NUMBER_OF_FILES << std::endl;
            }
            else
            {
                std::cout << "'" << argv[i] << "' is not a valid number." << std::endl;
            }
        }
        else if (std::strncmp(argv[i], "-alg",4) == 0)
        {
            i++;
            if (std::strncmp(argv[i], "isorank",7) == 0)
            {
                G_USE_ISORANK = true;
                std::cout << "Algorithm was set to IsoRank." << std::endl;
            }
            else if (std::strncmp(argv[i], "gpgm",4) == 0)
            {
                G_USE_GPGM = true;
                G_USE_ISORANK = false;
                std::cout << "Algorithm was set to GPGM." << std::endl;
            }
            else
            {
                std::cout << "Alorithm '" << argv [i] <<  "' is not a valid alorithm." << std::endl;
            }
        }
        else if (std::strncmp(argv[i], "-match-alg",7) == 0)
        {
            i++;
            if (std::strncmp(argv[i], "greedy",6) == 0)
            {
                G_USE_GREEDY_ALG = true;
                std::cout << "Matching algorithm was set to greedy." << std::endl;
            }
            else if (std::strncmp(argv[i], "con-enf-1", 9) == 0)
            {
                G_USE_CON_ENF_1 = true;
                G_USE_GREEDY_ALG = false;
                std::cout << "Matching algorithm was set to connectivity enforcement 1." << std::endl;
            }
            else if (std::strncmp(argv[i], "con-enf-2", 9) == 0)
            {
                G_USE_CON_ENF_2 = true;
                G_USE_GREEDY_ALG = false;
                std::cout << "Matching algorithm was set to connectivity enforcement 2." << std::endl;
            }
            else if (std::strncmp(argv[i], "con-enf-3", 9) == 0)
            {
                G_USE_CON_ENF_3 = true;
                G_USE_GREEDY_ALG = false;
                std::cout << "Matching algorithm was set to connectivity enforcement 3." << std::endl;
            }
            else if (std::strncmp(argv[i], "con-enf-4", 9) == 0)
            {
            	G_USE_CON_ENF_4 = true;
                G_USE_GREEDY_ALG = false;
                std::cout << "Matching algorithm was set to connectivity enforcement 4." << std::endl;
            }
            else
            {
                std::cout << "Alorithm '" << argv [i] <<  "' is not a valid alorithm." << std::endl;
            }
        }
        else
        {
            std::cout << "Pram '" << argv [i] <<  "' is not a valid argument." << std::endl;
        }
    }
    
    std::cout << "\n\n" <<"Program configuration: " << std::endl;
    std::cout << "Working directory: '" << G_DIR_PATH << "'" << std::endl;
    std::cout << "File Extention: '" << G_FILE_EXTENSION << "'" << std::endl;
    std::cout << "Number of graphs to read: " << G_NUMBER_OF_FILES << std::endl;
    if (G_USE_ISORANK)
    {
        std::cout << "Graph matching algorithm: IsoRank." << std::endl;
        std::string assignment_app;
        if (G_USE_GREEDY_ALG)
        {
            assignment_app = "Greedy.";
        }
        else if (G_USE_CON_ENF_1)
        {
            assignment_app = "Connectivity Enforcement 1.";
        }
        else if (G_USE_CON_ENF_2)
        {
            assignment_app = "Connectivity Enforcement 2.";
        }
        else if (G_USE_CON_ENF_3)
        {
            assignment_app = "Connectivity Enforcement 3.";
        }
        else if (G_USE_CON_ENF_4)
        {
            assignment_app = "Connectivity Enforcement 4.";
        }
        cout << "Assignment approach: "<< assignment_app << endl;
    }
    else if (G_USE_GPGM)
    {
        std::cout << "Graph matching algorithm: GPGM." << std::endl;
    }
    
    cout << '\n' << endl;
}
