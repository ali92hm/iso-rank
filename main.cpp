
/***********************************************************************************************************************************************
 * This is a software that gives an approximate solution to the graph matching problem by using the                                            *
 * IsoRank Algorithm. Created in the summer of 2013 this software was originally written with the                                              *
 * intention of being used by the Ferguson group at the Materials Science & Engineering Department at UIUC.                                    *
 *                                                                                                                                             *
 *                                                                                                                                             *
 * For more details on the software (i.e. design decisions, potential bugs etc.) please consult the README file.                               *
 *                                                                                                                                             *
 * Contact for Questions:                                                                                                                      *
 * Abhijit Pujare: abhijitpujare@gmail.com                                                                                                     *
 * Ali Hajimirza:  ali92hm@gmail.com                                                                                                           *
 *                                                                                                                                             *
 *                                                                                                                                             *
 *                                                                                                                                             *
 * Refs: M. Leordeanu and M. Herbert, "A Spectral Technique for Correspondence Problems Using Pairwise Constraints"                            *
 * Proceedings of the Tenth IEEE International Conference on Computer Vision (ICCV?05) 1550                                                    *
 *                                                                                                                                             *
 * R. Singh, J. Xu and B. Berger "Global alignment of multiple protein interaction networks with application to functional orthology detection"*
 * PNAS 105 35 12763?12768                                                                                                                     *
 *                                                                                                                                             *
 *                                                                                                                                             *
 * Special thanks to Dr. Andrew Ferguson and Andrew Long for giving us advice and guidance while creating this software.                       *
 *                                                                                                                                             *
 ***********************************************************************************************************************************************/


#include "Matrices/SymMatrix.h"
#include "Matrices/DenseMatrix2D.h"
#include "IsoRank.h"
#include <string>
#include <sstream>
#include <vector>
#include <ctime>
#include "mpi.h"

/*
 * Path to the graph folders.
 * Extension of the files.
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
 * Preprocessor definitions for the used data type
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


/*
 * Main function
 * @pram int argc
 * @pram char** argv
 */
int main(int argc, char * argv[])
{
    srand(time(NULL));
	/*
	 * MPI Variables
	 */
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
		total_comparisons = (0.5*(input_graphs.size()-1)*input_graphs.size());
		time_end = std::clock();
		if(G_PRINT)
			std::cout << input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
			<< timeElapsed(time_start, time_end) << "(ms)." << std::endl;
        
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
					input_graphs[i]->MPI_Send_SymMatrix(dest_ID, TAG_1 * dest_ID);
					input_graphs[j]->MPI_Send_SymMatrix(dest_ID, TAG_1 * dest_ID + TAG_2);
					
					if(G_DEBUG)
						std::cout <<"Master: sending matrix to ID: " << dest_ID << std::endl;
                    
					dest_ID++;
				}
				// Send additional pairs upon worker node's request
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
					input_graphs[i]->MPI_Send_SymMatrix(dest, TAG_1 * dest);
					input_graphs[j]->MPI_Send_SymMatrix(dest, TAG_1 * dest + TAG_2);
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
 			SymMatrix<DataType> emptyMat(0);
 			emptyMat.MPI_Send_SymMatrix (i, TAG_1*i);
 			emptyMat.MPI_Send_SymMatrix (i, TAG_1*i+ TAG_2);
 			if (G_DEBUG)
 				std::cout <<"Master: sending terminate signal to ID: " << i << std::endl;
 		}
 		time_end = std::clock();
        
		//printing the results
		if (G_PRINT)
		{
			std::cout << "Master: Computed IsoRank successfully for " << G_NUMBER_OF_FILES << " graphs in "
            << timeElapsed(time_start, time_end) << "(ms)." << std::endl;
            
			std::cout<< "Master: " <<isoRank_results.size() << " results were received.\n frob_norms: ";
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
    	while(true)
    	{
    		//Recv graphs from the master
    		DenseMatrix2D<DataType> mat1 (MASTER_ID, TAG_1 * ID ,stat);
    		DenseMatrix2D<DataType> mat2 (MASTER_ID, TAG_1 * ID + TAG_2 ,stat);
            
    		if (G_DEBUG)
    			std::cout << "Process "<< ID << " : received graphs from master"<< std::endl;
            
    		//Terminating the while loop if the matrices are 0*0
			if (mat1.getNumberOfRows() == 0 && mat2.getNumberOfRows() == 0)
			{
				if (G_DEBUG)
					std::cout << "Process "<< ID << ": received terminate signal from master"<< std::endl;
				break;
			}

			struct IsoRank_Result result;
			try
			{
				if (G_USE_ISORANK)
				{
					if (G_DEBUG)
						std::cout << "Process " << ID << ": isoRank: started." << std::endl;
			  		result = isoRank(mat1, mat2, G_GRAPH_MATCHING_ALGORITHM);
			  		if (G_DEBUG)
						std::cout << "Process " << ID << ": isoRank: end." << std::endl;
				}
				if (G_USE_GPGM)
				{
					//GPGM(mat1,mat2);
				}
                
				if (G_DEBUG)
					std::cout << "Process "<< ID << " :requesting for more graphs from master" << std::endl;
				//Sending the ID to master for more graphs
				MPI_Send(&ID, 1, MPI_INT, MASTER_ID, TAG_1 + TAG_2, MPI_COMM_WORLD);
                
				if (G_DEBUG)
                    std::cout << "Process "<< ID << ": :sending result to master" << std::endl;
				//Sending results to master
				MPI_Send_IsoRank_Result(result, MASTER_ID , TAG_1 * ID + TAG_3);
                delete [] result.assignments;
				
			}
			catch (std::exception& e)
			{
				std::cerr << "Process "<< ID << " Exception: " << e.what() << std::endl;
			}
		}
		
	}
	
	if(G_PRINT)
    	std::cout << "Process: "<< ID << " terminated." << std::endl;
    
    MPI_Finalize();
    return 0;
}
