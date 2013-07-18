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


int G_NUMBER_OF_FILES = 25;
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




void parseCommandLineArgs(int argc, char * argv[], int rank);
template <typename DT>
void MPI_Send_Matrix (SparseMatrix<DT>& matrix, int dest, int tag);
SparseMatrix<DataType> MPI_Recv_Matrix (int source, int tag , MPI_Status& stat);
//void  MPI_Recv_Matrix (int source, int tag , MPI_Status& stat);
struct IsoRank_Result* MPI_Recv_IsoRank_Result(int source, int tag, MPI_Status& stat);
void MPI_Send_IsoRank_Result (struct IsoRank_Result* result, int dest, int tag);



/*
 *
 */
int main(int argc, char * argv[])
{
  srand(time(NULL));
  
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
    parseCommandLineArgs(argc, argv, rank);
    std::vector<SparseMatrix<DataType> > input_graphs;
    std::vector<struct IsoRank_Results*> isoRank_results;
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
    	std::ostringstream itos_converter;
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

		time_start = std::clock();
		int counter = 2;
		for (int i = 0; i < input_graphs.size(); i++)
		{
			for(int j = i +1; j <  input_graphs.size(); j++)
			{
				if (counter < num_procs - 1)
				{
					MPI_Send_Matrix (input_graphs[i], counter, tag1*counter);
					MPI_Send_Matrix (input_graphs[j], counter, tag1*counter+10);
// 					std::cout <<"Master: sending matrix to rank:" << counter << std::endl;
					counter++;
				}
				else
				{
					int dest;
					MPI_Recv(&dest, 1, MPI_INT, MPI_ANY_SOURCE, tag1+10, MPI_COMM_WORLD,&stat);
					isoRank_results.push_back(MPI_Recv_IsoRank_Result(dest, tag1+11, &stat));
// 					std::cout <<"Master: received request for more graphs from: "<< dest<< std::endl;
					MPI_Send_Matrix (input_graphs[i], dest, tag1*dest);
					MPI_Send_Matrix (input_graphs[j], dest, tag1*dest+10);
// 					std::cout <<"Master: sending more graphs to: "<< dest<< std::endl;
				}
			}
		}
		
 		for(int i=1; i < num_procs -1; i++)
 		{
 			SparseMatrix<DataType> emptyMat(0,0);
 			MPI_Send_Matrix (emptyMat, i, tag1*i);
 			MPI_Send_Matrix (emptyMat, i, tag1*i+10);
//  			std::cout <<"Master: sending terminate signal to rank:" << i << std::endl;
 			
 		} 

    }
    else if (rank == 1)
    {
    	//collecting results
    }
    else
    {
    	while(true)
    	{
    		SparseMatrix<DataType> mat1 = MPI_Recv_Matrix (0, tag1*rank ,stat);
    		SparseMatrix<DataType> mat2 = MPI_Recv_Matrix (0, tag1*rank +10 ,stat);
			if (mat1.getNumberOfRows()==0 && mat2.getNumberOfRows()==0)
			{
// 				std::cout << "Process: "<< rank << "received terminate signal from master"<< std::endl;
				break;
			}
			else if(mat1.getNumberOfRows()<mat2.getNumberOfRows()){
			SparseMatrix<DataType>* mat_hold=&mat1;
			mat1=mat2;
			mat2=*mat_hold;
			
			}
			
		
			/*
			 * Compute n choose 2 combination of graphs
			 */

			
				int* assignment = new int[mat1.getNumberOfRows()];
				init_array(assignment,mat1.getNumberOfRows(),-1);
				struct IsoRank_Result* result;
					try
					{
						if (G_USE_ISORANK)
						{
							if (G_USE_GREEDY_ALG)
							{
								
					  			isoRank(mat1, mat2, 0,assignment);
// 								std::cout << "Process " << rank << ": isoRank " << endl;
							}
					
							else if (G_USE_CON_ENF_1)
							{
					 		 	isoRank(mat1, mat2, 1,assignment);
							}
					
							else if (G_USE_CON_ENF_2)
							{
					  			isoRank(mat1, mat2, 2,assignment);
							}
					
							else if (G_USE_CON_ENF_3)
							{
					  			isoRank(mat1, mat2, 3,assignment);
							}
							else if (G_USE_CON_ENF_4)
							{
					 			 isoRank(mat1, mat2, 4,assignment);
							}
						}
				
						if (G_USE_GPGM)
						{
							//GPGM(mat1,mat2);
						}

		
					}
					catch (std::exception& e)
					{
						std::cerr << "Exception: " << e.what() << std::endl;
					}
						

// 			  std::cout << "Process: "<< rank << ": requesting for more graphs from master" << std::endl;
			  
			  MPI_Send(&rank, 1, MPI_INT, 0, tag1+10, MPI_COMM_WORLD);
			  MPI_Recv_IsoRank_Result(result,0 , tag1+11);
			}
	}
	
	

//     std::cout << "Process: "<< rank << " terminated" << std::endl; 
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
    delete [] dataArray;
    return matrix;
}

void MPI_Send_IsoRank_Result (struct IsoRank_Result* result, int dest, int tag)
{
	MPI_Send(&rerult->length, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
	MPI_Send(&rerult->assignments, rerult->length, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
	MPI_Send(&rerult->forb_norm, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
}

struct IsoRank_Result* MPI_Recv_IsoRank_Result(int source, int tag, MPI_Status& stat)
{
	IsoRank_Result* result = new IsoRank_Result;
	MPI_Recv(&result->length, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
	results->assignemts = new int[result->length];
	MPI_Recv(results->assignemts ,result->length , MPI_INT, source, tag + 3, MPI_COMM_WORLD, &stat);
	MPI_Recv(&result->forb_norm, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
	return results;
}

/*
 *
 */
void parseCommandLineArgs(int argc,char* argv[], int rank)
{
  //  std::cout << "Reading command line arguments..." << std::endl;
    
    if ( argc < 2)
    {
    //    std::cout <<"No arguments were passed." << std::endl;
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
    
    if (rank == 0)
    {
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
}
