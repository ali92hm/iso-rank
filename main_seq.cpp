
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

#include "Matrices/DenseMatrix1D.h"
#include "IsoRank.h"
#include <vector>
#include <ctime>

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
    /*
     *Configure the program to use the command line args
     */
    srand(time(NULL));
    parseCommandLineArgs(argc, argv, 0);
    
    /*
     * Timing Variables
     */
    std::clock_t time_start;
    std::clock_t time_end;
	
	/*
	 *	Input/Result containers
	 */
	int total_comparisons;
    std::vector<IsoRank_Result> isoRank_results;
    std::vector<DenseMatrix1D<DataType>* >input_graphs;
    
    if(G_PRINT)
        std::cout << "Reading " << G_NUMBER_OF_FILES << " graphs from: " << G_DIR_PATH << std::endl;
    std::ostringstream itos_converter;
    time_start = std::clock();
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
    total_comparisons = (0.5*(input_graphs.size()-1)*input_graphs.size());
    time_end = std::clock();
    if(G_PRINT)
        std::cout << input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "
        << timeElapsed(time_start, time_end) << "(ms)." << std::endl;
    
    time_start = std::clock();
    for (int i = 0; i < input_graphs.size(); i++)
    {
        for (int j = i +1; j < input_graphs.size(); j++)
        {
            try
            {
                if (G_USE_ISORANK)
                {
                    isoRank_results.push_back(isoRank(*input_graphs[i], *input_graphs[j], G_GRAPH_MATCHING_ALGORITHM));
                }
                if (G_USE_GPGM)
                {
                    //GPGM(mat1,mat2);
                }
                
            }
            catch (std::exception& e)
            {
                std::cerr << " Exception: " << e.what() << std::endl;
            }
        }
    }
    time_end = std::clock();
    
	
    //printing the results
    if (G_PRINT)
    {
        std::cout << "Computed IsoRank successfully for " << G_NUMBER_OF_FILES << " graphs in "
        << timeElapsed(time_start, time_end) << "(ms)." << std::endl;
        
        std::cout << "Frob_norms: ";
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
    
    typename std::vector<DenseMatrix1D<DataType>* >::iterator graph_it;
    for ( graph_it = input_graphs.begin() ; graph_it < input_graphs.end(); ++graph_it )
    {
        delete  *graph_it;
    }
    return 0;
}

