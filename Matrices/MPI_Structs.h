


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

/*
 * function used to send the IsoRank_Result struct between two processors
 * @pram: result struct of isorank
 * @pram: the destination processor
 * @pram: the tag used for the MPI calls
 */
void MPI_Send_IsoRank_Result (IsoRank_Result result, int dest, int tag)
{
    MPI_Send(&result.assignment_length, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
    MPI_Send(result.assignments, result.assignment_length, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
    MPI_Send(&result.frob_norm, 1, MPI_INT, dest, tag + 3, MPI_COMM_WORLD);
}

/*
 * function used to receive the IsoRank_Result struct
 * @pram: the source processor where this struct came from
 * @pram: the tag used by the MPI calls
 * @pram: the MPI_Status object used by the MPI_calls
 */
struct IsoRank_Result MPI_Recv_IsoRank_Result(int source, int tag, MPI_Status& stat)
{
    struct IsoRank_Result result;
    MPI_Recv(&result.assignment_length, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
    result.assignments = new int[result.assignment_length];
    MPI_Recv(result.assignments ,result.assignment_length , MPI_INT, source, tag + 2, MPI_COMM_WORLD, &stat);
    MPI_Recv(&result.frob_norm, 1, MPI_INT, source, tag + 3, MPI_COMM_WORLD, &stat);
    return result;
}

#endif
#endif