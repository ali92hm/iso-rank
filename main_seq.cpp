/***********************************************************************************************************************************************
 * This is a software that gives an approximate solution to the graph matching problem by using the                                            *
 * IsoRank Algorithm. Created in the summer of 2013 this software was originally written with the                                              *
 * intention of being used by the Ferguson group at the Materials Science & Engineering Department at UIUC.                                    *
 * This is a sequential implementation of the code. A parallel version can be found in main.cpp.                                               *
 *                                                                                                                                             *
 * For more details on the software (i.e. how to run the program,design decisions, potential bugs etc.) please consult the README file.        *
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

#include "Matricies/DenseMatrix.h"
#include "IsoRank.h"
#include <vector>
#include <ctime>
#include "Tarjan.h"
#include "vertex.h"
#include "greedy_algorithms.h"

#ifdef __linux__
std::string G_DIR_PATH = "/home/ali/lib/input/";
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




void parseCommandLineArgs(int argc, char * argv[], int rank);

/*
 * main method that reads in graphs and calls the IsoRank method
 * @param: number of command line arguments
 * @param: array of strings of command line arguments
 */
int main(int argc, char * argv[])
{
 	srand(time(NULL));
 	
    
 	parseCommandLineArgs(argc, argv, 0);
 	std::vector<DenseMatrix<DataType> > input_graphs;
    //   std::vector<IsoRank_Result> isoRank_results;
 	std::clock_t time_start;
 	std::clock_t time_end;
 	double elapsed_time;
    
    //read graphs from files and enter them into an array of DenseMatrix objects
 	std::cout << "Reading " << G_NUMBER_OF_FILES << " graphs from: " << G_DIR_PATH << endl;
 	time_start = std::clock();
 	std::ostringstream itos_converter;
 	for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
 	{
 		try
 		{
 			itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
 			input_graphs.push_back(DenseMatrix<DataType>(itos_converter.str()));
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
    
    //run isorank on all pairs of graphs that are read in
 	for (int i = 0; i < input_graphs.size(); i++)
 	{
 		DenseMatrix<DataType> mat_A =  input_graphs[i];
 		for (int j = i+1; j < input_graphs.size(); j++)
 		{
 			
 			DenseMatrix<DataType> mat_B = input_graphs[j];
 			
 			struct IsoRank_Result result;
 			try
 			{
 				if (G_USE_ISORANK)
 				{
 					if (G_USE_GREEDY_ALG)
 					{
 						
 						result = isoRank(mat_A,mat_B , 0);
 						std::cout << result.frob_norm << std::endl;
 						delete []result.assignments;
 					}
 					
 					else if (G_USE_CON_ENF_1)
 					{
 						isoRank(mat_A, mat_B, 1);
 					}
 					
 					else if (G_USE_CON_ENF_2)
 					{
 						isoRank(mat_A, mat_B, 2);
 					}
 					
 					else if (G_USE_CON_ENF_3)
 					{
 						isoRank(mat_A, mat_B, 3);
 					}
 					else if (G_USE_CON_ENF_4)
 					{
 						isoRank(mat_A, mat_B, 4);
 					}
 				}
 				
 				if (G_USE_GPGM)
 				{
                    //GPGM(input_graphs[i],input_graphs[j]);
 				}
 			}
 			catch (std::exception& e)
 			{
 				std::cerr << "Exception: " << e.what() << std::endl;
 				
 			}
            //delete []assignment;
 		}
 	}
 	
 	time_end = std::clock();
 	elapsed_time = (double) (time_end - time_start) / CLOCKS_PER_SEC * 1000.0;
 	std::cout << "Computed IsoRank successfully for " << G_NUMBER_OF_FILES << " graphs in "<< elapsed_time << "(ms)." << endl;
    
 	return 0;
}


/*
 *function that parses command line arguments
 *@param: number of command line arguments passed
 *@param: command line arguments as an array of strings
 *@param: rank of the process calling parseCommandLineArgs, only rank 0 should print to console
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
 				std::cout << "Algorithm '" << argv [i] <<  "' is not a valid alorithm." << std::endl;
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
 				std::cout << "Algorithm '" << argv [i] <<  "' is not a valid alorithm." << std::endl;
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
