//
//  main.cpp
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/11/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#include "SparseMatrix.h"
#include "IsoRank.h"
#include <vector>
#include <ctime>
#include "Tarjan.h"
#include "vertex.h"
#include "greedy_algorithms.h"




/*
 *
 */
#ifdef __linux__
std::string G_DIR_PATH = "/export/home/reu_share/input/";
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




void parseCommandLineArgs(int argc, const char * argv[]);

/*
 *
 */
int main(int argc, const char * argv[])
{
    /*
     *Configure the program to use the command line args
     */
    parseCommandLineArgs(argc, argv);
    
    /*
     * Struecutes that contain the input
     */
    std::cout << "Reading " << G_NUMBER_OF_FILES << " graphs from: " << G_DIR_PATH << endl;
    std::clock_t time_start = std::clock();
    std::vector<SparseMatrix<DataType>*> input_graphs;
    ostringstream itos_converter;
    for(int i = 1; i <= G_NUMBER_OF_FILES; i++)
    {
        try
        {
            itos_converter << G_DIR_PATH << i << G_FILE_EXTENSION;
            input_graphs.push_back(new SparseMatrix<DataType>(itos_converter.str()));
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
    
    std::clock_t time_end = std::clock();
	double elapsed_time = (double) (time_end - time_start) / CLOCKS_PER_SEC * 1000.0;
    std::cout << input_graphs.size() << " of " << G_NUMBER_OF_FILES << " graphs were successfully read in "<< elapsed_time << "(ms)." << endl;
    
    
    
    /*
     * Compute n choose 2 combination of graphs
     */
    for (int i = 0; i < input_graphs.size(); i++)
    {
        for(int j = (i + 1); j <  input_graphs.size(); j++)
        {
            try
            {
                if (G_USE_ISORANK)
                {
                    if (G_USE_GREEDY_ALG)
                    {
                        isoRank(*input_graphs[i], *input_graphs[j], 0);
                    }
                    
                    else if (G_USE_CON_ENF_1)
                    {
                        isoRank(*input_graphs[i], *input_graphs[j], 1);
                    }
                    
                    else if (G_USE_CON_ENF_2)
                    {
                        isoRank(*input_graphs[i], *input_graphs[j], 2);
                    }
                    
                    else if (G_USE_CON_ENF_3)
                    {
                        isoRank(*input_graphs[i], *input_graphs[j], 3);
                    }
                    else if (G_USE_CON_ENF_4)
                    {
                        isoRank(*input_graphs[i], *input_graphs[j], 4);
                    }
                }
                
                if (G_USE_GPGM)
                {
                    //GPGM(*input_graphs[i], *input_graphs[j]);
                }
            }
            catch (std::exception& e)
            {
                std::cerr << "Exception: " << e.what() << std::endl;
            }
        }
    }
    
    /*
     * Deleting objects
     */
    vector<SparseMatrix<DataType>*>::iterator it;
    for ( it = input_graphs.begin(); it < input_graphs.end(); ++it )
    {
        delete (*it);
    }
    
    return 0;
}

/*
 *
 */
void parseCommandLineArgs(int argc,const char* argv[])
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
                std::cout << "Matching algorithm was set to enforcement 2." << std::endl;
            }
            else if (std::strncmp(argv[i], "con-enf-3", 9) == 0)
            {
                G_USE_CON_ENF_3 = true;
                G_USE_GREEDY_ALG = false;
                std::cout << "Matching algorithm was set to enforcement 3." << std::endl;
            }
            else if (std::strncmp(argv[i], "", 9) == 0) //TODO: andrew's algorithm
            {
                G_USE_GREEDY_ALG = false;
                std::cout << "Matching algorithm was set to ." << std::endl;
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
        else if (G_USE_CON_ENF_3)   //TODO Andrew's alg
        {
            assignment_app = "Connectivity Enforcement 3.";
        }
        cout << "Assignment approach: "<< assignment_app << endl;
    }
    else if (G_USE_GPGM)
    {
        std::cout << "Graph matching algorithm: GPGM." << std::endl;
    }
    
    cout << '\n' << endl;
}
