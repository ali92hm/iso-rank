//
//  main.cpp
//  Graph_Matching
//
//  Created by Ali Hajimirza on 6/11/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#include "Matricies/SymMatrix.h"
#include "Matricies/SymSparseMatrix.h"
#include "IsoRank.h"
#include <vector>
#include <ctime>
#include "mpi.h"

/*
 * Path to the graph folders.
 */
#ifdef __linux__
std::string G_DIR_PATH = "/home/ali/workspace/ex/input/";
#elif defined __APPLE__
std::string G_DIR_PATH = "/Users/AliHM/Documents/Course Material/Summer 13 REU/graphs/";
#endif

std::string G_FILE_EXTENSION = ".dat";
int G_NUMBER_OF_FILES = 4;

/*
 * 
 */
bool G_USE_ISORANK = true;
bool G_USE_GPGM = false;
int G_GRAPH_MATCHING_ALGORITHM = 0;
bool G_PRINT = false;
bool G_DEBUG = false;


/*
 * Preprocess definitions for the used data type
 */
#if INT_MATRIX
typedef int DataType;
#elif DOUBLE_MATRIX
typedef double DataType;
#elif LONG_DOUBLE
typedef long double DataType;
#else
typedef float DataType;
#endif

/*
 * Function protorypes
 */
void parseCommandLineArgs(int argc, char * argv[], int ID);
double timeElapsed(std::clock_t start, std::clock_t end);


/*
 * Main function
 * @pram int argc
 * @pram char** argv
 */
int main(int argc, char * argv[])
{
	/*
	 * MPI Variables
	 */  
	int num_procs;
    int ID;
 	MPI_Status stat;

 	/*
 	 *	MPI constant Tags
 	 */
 	const int MASTER_ID = 0;
    const int TAG_1 = 4;
    const int TAG_2 = 10;
    const int TAG_3 = 15;

    /*
     * MPI Initilalization calls 
     */
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        std::cout << "Failed To Initialize MPI" << std::endl;
        //MPI_Abort();
    }
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        
    /*
     *Configure the program to use the command line args
     */
    parseCommandLineArgs(argc, argv, rank);
    
    /*
     * Timing Varialbles
     */
    std::clock_t time_start;
    std::clock_t time_end;
	
	/*
	 *	Result variables
	 */
	int total_comparisons; 
    std::vector<IsoRank_Result> isoRank_results;

//======================================================================*MASTER NODE*==============================================================================
    if (ID == MASTER_ID)
    {
    	if(G_PRINT)
    		std::cout << "Reading " << G_NUMBER_OF_FILES << " graphs from: " << G_DIR_PATH << endl;
    	time_start = std::clock();
    	std::ostringstream itos_converter;
    	std::vector<SymMatrix<DataType>* >input_graphs;

    	/*
    	 * Reading the grphs and storing them
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
		total_comparisons = (0.5*(input_graphs.size()-1)*input_graphs.size());
		time_end = std::clock();
		if(G_PRINT)
			std::cout << input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "<< timeElapsed(time_start, time_end) << "(ms)." << endl;

		/*
    	 * Sending the graphs to worker nodes.
    	 */
		time_start = std::clock();
		int dest_ID = 1;
		int recv_counter = 0;
		for (int i = 0; i < input_graphs.size(); i++)
		{
			for(int j = i + 1; j <  input_graphs.size(); j++)
			{
				// Send a pair of graphs to all the worker nodes
				if (dest_ID < num_procs)
				{
					input_graphs[i]->MPI_Send_Matrix(dest_ID, TAG_1 * dest_ID);
					input_graphs[j]->MPI_Send_Matrix(dest_ID, TAG_1 * dest_ID + TAG_2);
					dest_ID++;
					if(G_DEBUG)
						std::cout <<"Master: sending matrix to rank:" << dest_ID << std::endl;
				}
				// Send additional pais upon worker node's request
				else
				{
					int dest;
					//Recv workers ID
					MPI_Recv(&dest, 1, MPI_INT, MPI_ANY_SOURCE, TAG_1 + TAG_2, MPI_COMM_WORLD,&stat);						
					if(G_DEBUG)
						std::cout <<"Master: received request for more graphs from: "<< dest<< std::endl;

					//Collect the result from worker
					isoRank_results.push_back(MPI_Recv_IsoRank_Result(dest, TAG_1 * dest + TAG_3, stat));						
					recv_counter++;
					if(G_DEBUG)
						std::cout <<"Master: results were received "<< dest<< std::endl;

					//Send more graphs to worker node
					MPI_Send_Matrix (input_graphs[i], dest, TAG_1 * dest);
					MPI_Send_Matrix (input_graphs[j], dest, TAG_1 * dest + TAG_2);
					if(G_DEBUG)
						std::cout <<"Master: sending more graphs to: "<< dest<< std::endl;
				}
			}
		}
		// Recv the remaining result
		while (recv_counter < total_comparisons)
		{
			int dest;
			//Recv workers ID
			MPI_Recv(&dest, 1, MPI_INT, MPI_ANY_SOURCE, TAG_1 + TAG_2, MPI_COMM_WORLD,&stat);
			if(G_DEBUG)
				std::cout <<"Master: received request for more graphs from: "<< dest<< std::endl;

			//Collect the result from worker
			isoRank_results.push_back(MPI_Recv_IsoRank_Result(dest, TAG_1 * dest + TAG_3, stat));
			recv_counter++;
			if(G_DEBUG)
				std::cout <<"Master: results were received."<< dest<< std::endl;
		}
		
		//Terminating the slaves by sending a 0*0 matrix to nodes
 		for(int i=1; i < num_procs; i++)
 		{
 			SparseMatrix<DataType> emptyMat(0,0);
 			MPI_Send_Matrix (emptyMat, i, TAG_1*i);
 			MPI_Send_Matrix (emptyMat, i, TAG_1*i+10);
 			if (G_DEBUG)
 				std::cout <<"Master: sending terminate signal to rank:" << i << std::endl;	
 		} 
 		time_end = std::clock();

		if (G_PRINT)
			std::cout << "Computed IsoRank successfully for " << G_NUMBER_OF_FILES << " graphs in "<< timeElapsed(time_start, time_end) << "(ms)." << endl;

		//printing the results 
		if (G_PRINT)
		{
			for (int i=0; i < isoRank_results.size(); i++)
			{
				std::cout<< isoRank_results[i].frob_norm << ", ";
			}
		}
    }
//======================================================================*WORKER NODES*==============================================================================
    else
    {
    	while(true)
    	{
    		//Recv graphs from the master
    		SymSparseMatrix<DataType> mat1 = SymSparseMatrix (MASTER_ID, TAG_1 * ID ,stat);
    		SymSparseMatrix<DataType> mat2 = SymSparseMatrix (MASTER_ID, TAG_1 * ID + TAG_2 ,stat);

    		if (G_DEBUG)
    			std::cout << "Process: "<< rank << "received graphs from master"<< std::endl;

    		//Terminating the while loop if the matrices are 0*0
			if (mat1.getNumberOfRows() == 0 && mat2.getNumberOfRows() == 0)
			{
				if (G_DEBUG)
					std::cout << "Process: "<< rank << "received terminate signal from master"<< std::endl;
				break;
			}	

			int* assignment = new int[mat1.getNumberOfRows()];
			init_array(assignment,mat1.getNumberOfRows(),-1);
			try
			{
				if (G_USE_ISORANK)
				{	
					if (G_DEBUG)
						std::cout << "Process " << rank << ": isoRank: started." << endl;
			  			result = isoRank(mat1, mat2, 0,assignment);
			  		if (G_DEBUG)
						std::cout << "Process " << rank << ": isoRank: end." << endl;
				if (G_USE_GPGM)
				{
					//GPGM(mat1,mat2);
				}

				if (G_DEBUG)
					std::cout << "Process: "<< rank << ": requesting for more graphs from master" << std::endl;
				//Sending the ID to master for more graphs
				MPI_Send(&rank, 1, MPI_INT, MASTER_ID, TAG_1 + TAG_2, MPI_COMM_WORLD);

				if (G_DEBUG)
				  std::cout << "Process: "<< rank << ": sending result to master" << std::endl;
				//Sending results to master
				MPI_Send_IsoRank_Result(result, MASTER_ID , TAG_1 * ID + TAG_3);
				
			}
			catch (std::exception& e)
			{
				std::cerr << "Exception: " << e.what() << std::endl;
			}
		}			
		
	}
	
	if(G_PRINT)
    	std::cout << "Process: "<< ID << " terminated." << std::endl; 

    MPI_Finalize();
    return 0;
}

/*
 * Calculates the time elapsed
 * @pram std::clock_t  start_time
 * @pram std::clock_t  end_time
 */
double timeElapsed(std::clock_t start, std::clock_t end)
{
	return (double) (time_end - time_start) / CLOCKS_PER_SEC * 1000.0;
}

/*
 * This method Configures the setting of the program according to the command line agrs.
 * @pram int argc
 * @pram char** argv
 * @pram int processor ID
 */
void parseCommandLineArgs(int argc,char* argv[], int ID)
{
   if (ID == 0)
   {
	   std::cout << "Reading command line arguments..." << std::endl; 
	   if ( argc < 2)
	   {
	      std::cout <<"No arguments were passed." << std::endl;
	   }
	}
    
    for (int i = 1; i < argc; i++)
    {
    	//changing the folder directory
        if (std::strncmp(argv[i], "-dir",4) == 0)
        {
            i++;
            G_DIR_PATH = std::string(argv[i]);
            if (ID == 0)
            	std::cout << "Working directory was set to: " << G_DIR_PATH << std::endl;
        }
        //changing the extention of the file
        else if (std::strncmp(argv[i], "-ext",4) == 0)
        {
            i++;
            G_FILE_EXTENSION = std::string(argv[i]);
            if (ID == 0)
            	std::cout << "File extention was set to: " << G_FILE_EXTENSION << std::endl;
        }
        //changing the number of files to be read
        else if (std::strncmp(argv[i], "-num_files",10) == 0)
        {
            i++;
            int input_number = atoi(argv[i]);
            if ( input_number > 0)
            {
                G_NUMBER_OF_FILES = input_number;
                if (ID == 0)
                	std::cout << "Number of files was set to: " << G_NUMBER_OF_FILES << std::endl;
            }
            // the input is not a number or it's an invalid number
            else
            {
            	if (ID == 0)
                	std::cout << "'" << argv[i] << "' is not a valid number." << std::endl;
            }
        }
        //changing the Graph Matching algorithm
        else if (std::strncmp(argv[i], "-alg",4) == 0)
        {
            i++;
            //choosing IsoRnak
            if (std::strncmp(argv[i], "isorank",7) == 0)
            {
                G_USE_ISORANK = true;
                if (ID == 0)
                	std::cout << "Algorithm was set to IsoRank." << std::endl;
            }
            //choosing GPGM
            else if (std::strncmp(argv[i], "gpgm",4) == 0)
            {
                G_USE_GPGM = true;
                G_USE_ISORANK = false;
                if (ID == 0)
                	std::cout << "Algorithm was set to GPGM." << std::endl;
            }
            //case of invalid pram for algorithm
            else
            {
            	if (ID == 0)
                	std::cout << "Alorithm '" << argv [i] <<  "' is not a valid alorithm." << std::endl;
            }
        }
        //chanigng the score matching algorithm
        else if (std::strncmp(argv[i], "-match_alg",7) == 0)
        {
            i++;
            //choosing greedy
            if (std::strncmp(argv[i], "greedy",6) == 0)
            {
                G_GRAPH_MATCHING_ALGORITHM = 0;
                if (ID == 0)
                	std::cout << "Matching algorithm was set to greedy." << std::endl;
            }
            //choosing connectivity enforcement 1
            else if (std::strncmp(argv[i], "con_enf_1", 9) == 0)
            {
                G_GRAPH_MATCHING_ALGORITHM = 1;
                G_USE_GREEDY_ALG = false;
                if (ID == 0)
                	std::cout << "Matching algorithm was set to connectivity enforcement 1." << std::endl;
            }
            //choosing connectivity enforcement 2
            else if (std::strncmp(argv[i], "con_enf_2", 9) == 0)
            {
                G_GRAPH_MATCHING_ALGORITHM = 2;
                if (ID == 0)
                	std::cout << "Matching algorithm was set to connectivity enforcement 2." << std::endl;
            }
            //choosing connectivity enforcement 3
            else if (std::strncmp(argv[i], "con_enf_3", 9) == 0)
            {
                G_GRAPH_MATCHING_ALGORITHM = 3;
                if (ID == 0)
                	std::cout << "Matching algorithm was set to connectivity enforcement 3." << std::endl;
            }
             //choosing connectivity enforcement 4
            else if (std::strncmp(argv[i], "con_enf_4", 9) == 0)
            {
            	G_GRAPH_MATCHING_ALGORITHM = 4;
                if (ID == 0)
                	std::cout << "Matching algorithm was set to connectivity enforcement 4." << std::endl;
            }
            // invalid pram for matching alorithm
            else
            {
            	if (ID == 0)
                	std::cout << "Alorithm '" << argv [i] <<  "' is not a valid alorithm." << std::endl;
            }
        }
        //Print to console
        else if (std::strncmp(argv[i], "-print", 9) == 0)
        {
        	i++;
        	G_PRINT = true;
        	if (std::strncmp(argv[i], "debug", 5) == 0)
            	G_DEBUG = true;
            if (ID == 0)
                	std::cout << "Printing to console was enabled.'"<< std::endl;
        }
        // invalid pram
        else
        {
        	if (ID == 0)
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
			switch (G_GRAPH_MATCHING_ALGORITHM)
			{
				case 0: 
					assignment_app = "Greedy.";
					break;
				case 1:
					assignment_app = "Connectivity Enforcement 1.";
					break;
				case 2:
					assignment_app = "Connectivity Enforcement 2.";
					break;
				case 3:
					assignment_app = "Connectivity Enforcement 3.";
					break;
				case 4:
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
