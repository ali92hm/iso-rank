/************************************************************************************
 * This is a file that is used to perform the IsoRank Algorithm.                    *
 * The Kroencker Product and the eigenvalue decomposition are done in this          *
 * file to get the scores matrix between nodal pairs. The ARPACK++ library          *
 * was used to compute the eigenvector decomposition. Greedy algorithms             *
 * to do the matchings are called in this file. Furthermore, functions used to      *
 * send and receive results of IsoRank between processors are defined in this file. *
 *                                                                                  *
 ************************************************************************************/

#ifndef Sparse_Matrix_IsoRank_h
#define Sparse_Matrix_IsoRank_h

#include "Matricies/DenseMatrix2D.h"
#include "Tarjan.h"
#include "util.h"
#include "greedy_algorithms.h"
#include <vector>

const int GREEDY = 0;
const int CON_ENF_1 = 1;
const int CON_ENF_2 = 2;
const int CON_ENF_3 = 3;
const int CON_ENF_4 = 4;

const int NUM_OF_ISORANK_IT = 20;


/*
 * struct used to store the result
 * of each isorank comparison done between
 * two graphs
 */
struct IsoRank_Result
{
	int frob_norm;
    int assignment_length;
	int* assignments;
};

/*
 * function used to perform the isorank algorithm
 * @param: adjacency matrix for graph1
 * @param: adjacency matrix for graph2
 * @param: the matching algorithm used to choose the best node to node mapping
 */
template <typename T>
struct IsoRank_Result isoRank(DenseMatrix2D<T>& matrix_A, DenseMatrix2D<T>& matrix_B, int matching_algorithm)
{
    //check to see both adjacency matrices are square and symmetric
    if (!matrix_A.isSquare() || !matrix_B.isSquare())
    {
        throw NotASquareMatrixException();
    }
    
    if (!matrix_A.isSymmetric() || !matrix_B.isSymmetric())
    {
        throw NotASymmetricMatrixException();
    }
    
    // Degree distribution statistics
    DenseMatrix2D<T> kron_prod = matrix_A.kron(matrix_B);
    std::vector<vertex*> vertices = graph_con_com(kron_prod);
    T** eigenValues = new T*[kron_prod.getNumberOfColumns()];
    vector<vector<int>*> comp_mask_values (kron_prod.getNumberOfColumns());
    
    //for each component find the eigenvector corresponding to the scores matrix
    for(int i=0; i < kron_prod.getNumberOfColumns(); i++ )
    {
        vector<int>* comp_mask = component_mask(vertices, i);
        comp_mask_values[i] = comp_mask;
        if (comp_mask == NULL)
        {
	  		eigenValues[i]=NULL;
            continue;
        }
        
        DenseMatrix2D<T> L = kron_prod.getScatteredSelection(*comp_mask,*comp_mask);
        
        std::vector<T> sum = L.getSumOfRows();
        std::vector<T> D_neg1(sum.size());
        std::vector<T> D_0pt5(sum.size());
        std::vector<T> D_neg0pt5(sum.size());
        
        for(int j=0; j < sum.size(); j++)
        {
            D_neg1[j] = 1.0/sum[j];
            D_0pt5[j] = sqrt(sum[j]);
            D_neg0pt5[j] = 1.0/D_0pt5[j];
        }
        
        DenseMatrix2D<T> lTimesD = L.diagonalVectorTimesMatrix(D_neg0pt5);
        DenseMatrix2D<T> Ms = lTimesD.matrixTimesDiagonalVector(D_neg0pt5);
        
        if(!Ms.isSymmetric())
        {
            throw NotASymmetricMatrixException();
        }
        
        T* eigenVec=  Ms.getTopEigenVector();
        double vecLength = 0;
        for (int j=0; j < L.getNumberOfRows(); j++)
        {
            eigenVec[j] *= D_0pt5[j];
            vecLength += pow(eigenVec[j],2);
        }
        
        vecLength = sqrt(vecLength);
        
        int coef = 1;
        if( eigenVec[0] < 0)
        {
            coef = -1;
        }
        
        for (int j=0; j < L.getNumberOfRows(); j++)
        {
            eigenVec[j] = coef* (eigenVec[j]/vecLength);
        }
        eigenValues[i] = eigenVec;
    }
    
    
    DenseMatrix2D<T> scores;
    vector<int>* comp_mask_curr=comp_mask_values[0];
    int counter_eig_vector=0;
    int counter_comp_mask=0;
    struct IsoRank_Result ret_val;
    
    //run matching algorithms for each component
    for(int k=0;k<kron_prod.getNumberOfColumns();k++) {
        T* eigenvector=eigenValues[k];
        
        if(eigenvector!=NULL) {
            comp_mask_curr=comp_mask_values[k];
            scores= reshape(eigenvector,matrix_A.getNumberOfRows(),matrix_B.getNumberOfColumns(),*comp_mask_curr);
            DenseMatrix2D<T> scores_copy(scores);
            int * best_assignment;
            float best_frob_norm=DBL_MAX;
            
            for (int num_it = 0; num_it < NUM_OF_ISORANK_IT; num_it++){
                scores=scores_copy;
                int* assignment= new int[matrix_A.getNumberOfRows()];
                init_array(assignment,matrix_A.getNumberOfRows(),-1);
                
                //call the approprite matching algorithm
                switch (matching_algorithm)
                {
                    case GREEDY:
                        greedy_1(scores,matrix_A,matrix_B,assignment);
                        break;
                    case CON_ENF_1:
                        greedy_connectivity_1(scores,matrix_A,matrix_B,assignment);
                        break;
                    case CON_ENF_2:
                        greedy_connectivity_2(scores,matrix_A,matrix_B,assignment);
                        break;
                    case CON_ENF_3:
                        greedy_connectivity_3(scores,matrix_A,matrix_B,assignment);
                        break;
                    case CON_ENF_4:
                        greedy_connectivity_4(scores,matrix_A,matrix_B,assignment);
                        break;
                    default:
                        break;
                }
                
                
                //find the frobenius norm
                // two cases: matrixA is larger than matrixB and matrixA is smaller than matrixB
                if(matrix_A.getNumberOfRows()>matrix_B.getNumberOfRows()){
                    DenseMatrix2D<T> perm_mat=getPermMatrix(assignment,matrix_A.getNumberOfRows(),matrix_A.getNumberOfRows());
                    DenseMatrix2D<T> product=perm_mat*matrix_A;
                    DenseMatrix2D<T> get_transpose=perm_mat.transpose();
                    DenseMatrix2D<T> final_mat=product*get_transpose;
                    DenseMatrix2D<T> ret_matrix= matrix_A-final_mat;
					
                    T frob_norm_hold=ret_matrix.getFrobNorm();
                    
                    if(frob_norm_hold<best_frob_norm){
                        best_frob_norm=frob_norm_hold;
                        if(num_it>0)
                            delete []best_assignment;
                        best_assignment=assignment;
                        
                    }
                    else{
                        delete []assignment;
                    }
                }
                else{
                    int b_size=matrix_B.getNumberOfRows();
                    int a_size=matrix_A.getNumberOfRows();
                    DenseMatrix2D<float> matrix_A2(b_size,b_size);
                    
                    for(int r=0;r<b_size;r++){
                        for(int s=0;s<b_size;s++){
                            if(r<a_size&&s<a_size){
                                matrix_A2(r,s)=matrix_A(r,s);
                            }
                            else if(r==s){
                                matrix_A2(r,s)=1;
                            }
                            else{
                                matrix_A2(r,s)=0;
                            }
                        }
                    }
                    
                    DenseMatrix2D<T> perm_mat=getPermMatrix(assignment,matrix_A.getNumberOfRows(),matrix_B.getNumberOfRows());
                    DenseMatrix2D<T> product=perm_mat*matrix_A2;	
                    DenseMatrix2D<T> get_transpose=perm_mat.transpose();
                    DenseMatrix2D<T> final_mat=product*get_transpose;
                    DenseMatrix2D<T> ret_matrix= matrix_A-final_mat;
                    
                    float frob_norm_hold=ret_matrix.getFrobNorm(); 
                    
                    if(frob_norm_hold<best_frob_norm){
                        best_frob_norm=frob_norm_hold;
                        if(num_it>0)
                            delete []best_assignment;
                        best_assignment=assignment;
                    }
                    else{
                        delete []assignment;
                    }
                }
            }
            
            ret_val.frob_norm=best_frob_norm;
            ret_val.assignments=best_assignment;
            ret_val.assignment_length=matrix_A.getNumberOfRows();        
        }
    }
    
	for(int k=0;k<kron_prod.getNumberOfColumns();k++)
	{
        delete []eigenValues[k];
	}
    
    for (int i = 0 ; i < comp_mask_values.size(); ++i )
    {
        delete  comp_mask_values[i];
    }
    
    for (int i=0; i < vertices.size() ; i++)
    {
    	delete vertices[i];
    }
    
    delete []eigenValues;
    return ret_val;
}

#ifdef USE_MPI
/*
 * function used to send the IsoRank_Result struct between two processors
 * @param: result struct of isorank
 * @param: the destination processor
 * @param: the tag used for the MPI calls
 */
void MPI_Send_IsoRank_Result (IsoRank_Result result, int dest, int tag)
{
	MPI_Send(&result.assignment_length, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
	MPI_Send(result.assignments, result.assignment_length, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
	MPI_Send(&result.frob_norm, 1, MPI_INT, dest, tag + 3, MPI_COMM_WORLD);
}

/*
 * function used to receive the IsoRank_Result struct
 * @param: the source processor where this struct came from
 * @param: the tag used by the MPI calls
 * @param: the MPI_Status object used by the MPI_calls
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
