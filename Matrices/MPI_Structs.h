


#ifndef __MPI_Structs__
#define __MPI_Structs__
#ifdef USE_MPI
#include "mpi.h"


static const unsigned short _DENSE_FORM = 0;
static const unsigned short _SPARSE_FORM = 1;
static const unsigned short _SYM_DENSE_FORM = 2;
static const unsigned short _SYM_SPARSE_FORM = 3;

struct MPI_MatrixInfo
{
    unsigned short rows;
    unsigned short cols;
    unsigned short recv_size;
    unsigned short send_form;
    void setValues(unsigned short num_rows, unsigned short num_cols, unsigned short num_recv_size, const unsigned short const_send_form)
    {
        rows = num_rows;
        cols = num_cols;
        recv_size = num_recv_size;
        send_form = const_send_form;
    }
};

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

#endif
#endif