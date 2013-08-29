//
//  main.cpp
//  Graph_Matching
//
//  Created by Ali Hajimirza on 6/11/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#include "Matrices/SymMatrix.h"
#include "Matrices/DenseMatrix2D.h"
#include "IsoRank.h"
#include <vector>
#include <ctime>
#include "mpi.h"

/*
 * Path to the graph folders.
 * Extention of the files.
 * Number of files to read.
 */
std::string G_DIR_PATH = "Sample input/";
std::string G_FILE_EXTENSION = ".dat";
int G_NUMBER_OF_FILES = 2;

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
 * Function prototypes
 */
void parseCommandLineArgs(int argc, char * argv[], int ID);
double timeElapsed(std::clock_t start, std::clock_t end);

struct Offset
{
    int i_start;
    int i_end;
    int j_start;
    int j_end;
    void setValues(int ID, int num_procs, int number_of_graphs)
    {
        int number_of_comparisons = number_of_graphs*(number_of_graphs-1)/2;
        int block_size = std::ceil((float)number_of_comparisons/(float)(num_procs-1));
        int block_start = (ID - 1) * block_size + 1;
        int block_end = std::min(ID*block_size, number_of_comparisons);
        i_start = std::floor(((2*number_of_graphs-1) - std::sqrt(9-4*number_of_graphs*(1-number_of_graphs)-8*block_start))/2);
        i_end = std::floor(((2*number_of_graphs-1) - std::sqrt(9-4*number_of_graphs*(1-number_of_graphs)-8*block_end))/2);
        j_start = block_start + (i_start*(i_start+3-2*number_of_graphs))/2;
        j_end = block_end + (i_end*(i_end+3-2*number_of_graphs))/2;
    }

};

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
     * MPI Initialization calls 
     */
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        std::cout << "Failed To Initialize MPI" << std::endl;
        //MPI_Abort();
    }
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank (MPI_COMM_WORLD, &ID);
        
    /*
     *Configure the program to use the command line args
     */
    parseCommandLineArgs(argc, argv, ID);
    
    /*
     * Timing Variables
     */
    std::clock_t time_start;
    std::clock_t time_end;
	
	/*
	 *	Result variables
	 */
	int total_comparisons; 
	int number_of_graphs;
    std::vector<IsoRank_Result> isoRank_results;

//======================================================================*MASTER NODE*==============================================================================
    if (ID == MASTER_ID)
    {
    	if(G_PRINT)
    		std::cout << "Reading " << G_NUMBER_OF_FILES << " graphs from: " << G_DIR_PATH << std::endl;
    	time_start = std::clock();
    	std::ostringstream itos_converter;
    	std::vector<SymMatrix<DataType>* >input_graphs;

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
		number_of_graphs = input_graphs.size();
		total_comparisons = (0.5*(number_of_graphs-1)*number_of_graphs);
		time_end = std::clock();
		if(G_PRINT)
			std::cout << input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
			<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;

		/*
    	 * Sending the graphs to worker nodes.
    	 */
		time_start = std::clock();
		MPI_Bcast (&number_of_graphs, 1 , MPI_INT, MASTER_ID, MPI_COMM_WORLD);
		if(G_DEBUG)
			std::cout <<"Master: sending "<<number_of_graphs << " graphs to all"<< std::endl;
		for (int i = 0; i < input_graphs.size(); i++)
		{
			input_graphs[i]->MPI_Bcast_Send_SymMatrix(MASTER_ID);

		}
		
		/*
		 * Collecting the results from the worker nodes.
		 */
		int recv_counter = 0;
		while (recv_counter < total_comparisons)
		{
			int dest;
			//Recv workers ID
			MPI_Recv(&dest, 1, MPI_INT, MPI_ANY_SOURCE, TAG_1 + TAG_2, MPI_COMM_WORLD,&stat);
			if(G_DEBUG)
				std::cout <<"Master: received signal to receive result from: "<< dest<< std::endl;

			//Collect the result from worker
			isoRank_results.push_back(MPI_Recv_IsoRank_Result(dest, TAG_1 * dest + TAG_3, stat));
			recv_counter++;
			if(G_DEBUG)
				std::cout <<"Master: results were received."<< dest<< std::endl;
		}
		
 		time_end = std::clock();

		//printing the results 
		if (G_PRINT)
		{
			std::cout << "Master: Computed IsoRank successfully for " << G_NUMBER_OF_FILES << " graphs in "
												<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;
												
			std::cout<< "Master: " << isoRank_results.size() << " results were received.\n frob_norms: ";
			for (int i=0; i < isoRank_results.size(); i++)
			{
				std::cout<< isoRank_results[i].frob_norm << ", ";
			}
			std::cout<<std::endl;
		}
		
		typename std::vector<IsoRank_Result>::iterator res_it;
		for ( res_it = isoRank_results.begin() ; res_it < isoRank_results.end(); ++res_it )
		{
			delete [] res_it->assignments;
		}
		
		typename std::vector<SymMatrix<DataType>* >::iterator graph_it;
		for ( graph_it = input_graphs.begin() ; graph_it < input_graphs.end(); ++graph_it )
		{
			delete  *graph_it;
		}
    }
//======================================================================*WORKER NODES*==============================================================================
    else
    {
    	
    	std::vector<DenseMatrix2D<DataType>* > recv_graphs;
    	
    	MPI_Bcast (&number_of_graphs, 1, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    	for (int i = 0; i < number_of_graphs; i++)
    	{
    		recv_graphs.push_back(new DenseMatrix2D<DataType>(MASTER_ID, stat));
    	}

		if (G_DEBUG)
			std::cout << "Process "<< ID << " : received " << number_of_graphs << " graphs from master"<< std::endl;
	
        Offset offset;
        offset.setValues(ID, num_procs, number_of_graphs);
        int A, B;
		for (int i = offset.i_start; i <= offset.i_end; i++)
		{	
            if (i == offset.i_start)
            {
                A = offset.j_start;
            }
            else
            {
                A = i+1;
            }

            if (i == offset.i_end)
            {
                B = offset.j_end;
            }
            else 
            {
                B = number_of_graphs - 1;
            }

			for (int j = A; j <= B; j++)
			{	
				struct IsoRank_Result result;
				try
				{
					if (G_USE_ISORANK)
					{	
						if (G_DEBUG)
							std::cout << "Process " << ID << ": isoRank: started."  << i << " " << j << std::endl;
						result = isoRank(*recv_graphs[i], *recv_graphs[j], G_GRAPH_MATCHING_ALGORITHM);
						if (G_DEBUG)
							std::cout << "Process " << ID << ": isoRank: end." << std::endl;
					}
					if (G_USE_GPGM)
					{
						//GPGM(mat1,mat2);
					}
						
					//Sending the ID to master for more graphs
					MPI_Send(&ID, 1, MPI_INT, MASTER_ID, TAG_1 + TAG_2, MPI_COMM_WORLD);

					if (G_DEBUG)
					  std::cout << "Process "<< ID << " :sending result to master" << std::endl;
					  
					//Sending results to master
					MPI_Send_IsoRank_Result(result, MASTER_ID , TAG_1 * ID + TAG_3);
					delete []result.assignments;
				}
				catch (std::exception& e)
				{
					std::cerr << "Process "<< ID << " Exception: " << e.what() << std::endl;
				}
			}
		}
		typename std::vector<DenseMatrix2D<DataType>* >::iterator graph_it;
		for ( graph_it = recv_graphs.begin() ; graph_it < recv_graphs.end(); ++graph_it )
		{
			delete  *graph_it;
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
	return (double) (end - start) / CLOCKS_PER_SEC * 1000.0;
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
        std::cout << "\n\nReading command line arguments..." << std::endl;
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
        //changing the extension of the file
        else if (std::strncmp(argv[i], "-ext",4) == 0)
        {
            i++;
            G_FILE_EXTENSION = std::string(argv[i]);
            if (ID == 0)
                std::cout << "File extension was set to: " << G_FILE_EXTENSION << std::endl;
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
                    std::cout << "Algorithm '" << argv [i] <<  "' is not a valid algorithm." << std::endl;
            }
        }
        //changing the score matching algorithm
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
            // invalid pram for matching algorithm
            else
            {
                if (ID == 0)
                    std::cout << "Algorithm '" << argv [i] <<  "' is not a valid algorithm." << std::endl;
            }
        }
        //Print to console
        else if (std::strncmp(argv[i], "-print", 6) == 0)
        {
            G_PRINT = true;
            if (ID == 0)
                std::cout << "Printing to console: enabled."<< std::endl;
        }
        //Print debug to console
        else if (std::strncmp(argv[i], "-debug", 6) == 0)
        {
            G_DEBUG = true;
            if (ID == 0)
                std::cout << "Debug prints: enabled."<< std::endl;
        }
        // invalid pram
        else
        {
            if (ID == 0)
                std::cout << "Arg '" << argv [i] <<  "' is not a valid argument." << std::endl;
        }
    }
    if (ID == 0)
    {
        std::cout << "\n\n" <<"Program configuration: " << std::endl;
        std::cout << "Working directory: '" << G_DIR_PATH << "'" << std::endl;
        std::cout << "File Extension: '" << G_FILE_EXTENSION << "'" << std::endl;
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
            std::cout << "Assignment approach: "<< assignment_app << std::endl;
        }
        else if (G_USE_GPGM)
        {
            std::cout << "Graph matching algorithm: GPGM." << std::endl;
        }
        std::cout << '\n' << std::endl;
    }
}