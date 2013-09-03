/************************************************************************************
 * This is a file that is used to perform the IsoRank Algorithm.                    *
 * The Kroencker Product and the eigenvalue decomposition are done in this          *
 * file to get the scores matrix between nodal pairs. The ARPACK++ library          *
 * was used to compute the eigenvector decomposition. Greedy algorithms             *
 * to do the matchings are called in this file. Furthermore, functions used to      *
 * send and receive results of IsoRank between processors are defined in this file. *
 *                                                                                  *
 ************************************************************************************/

#ifndef _IsoRank_h
#define _IsoRank_h

#include "Matrices/DenseMatrix1D.h"
#include "Matrices/MPI_Structs.h"
#include "Tarjan.h"
#include "Utilities.h"
#include "GreedyAlgorithms.h"
#include <vector>


static const int GREEDY = 0;
static const int CON_ENF_1 = 1;
static const int CON_ENF_2 = 2;
static const int CON_ENF_3 = 3;
static const int CON_ENF_4 = 4;

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
 * @pram: adjacency matrix for graph1
 * @pram: adjacency matrix for graph2
 * @pram: the matching algorithm used to choose the best node to node mapping
 */
template <typename T>
struct IsoRank_Result isoRank(DenseMatrix1D<T>& matrix_A, DenseMatrix1D<T>& matrix_B, int matching_algorithm)
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
    DenseMatrix1D<T> kron_prod = matrix_A.kron(matrix_B);
    std::vector<vertex*> vertices = graph_con_com(kron_prod);
    T** eigenValues = new T*[kron_prod.getNumberOfColumns()];
    std::vector<std::vector<int>*> comp_mask_values (kron_prod.getNumberOfColumns());
    
    //for each component find the eigenvector corresponding to the scores matrix
    for(int i=0; i < kron_prod.getNumberOfColumns(); i++ )
    {
        std::vector<int>* comp_mask = component_mask(vertices, i);
        comp_mask_values[i] = comp_mask;
        if (comp_mask == NULL)
        {
	  		eigenValues[i]=NULL;
            continue;
        }
        
        DenseMatrix1D<T> L = kron_prod.getScatteredSelection(*comp_mask,*comp_mask);
        
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
        
        DenseMatrix1D<T> lTimesD = L.diagonalVectorTimesMatrix(D_neg0pt5);
        DenseMatrix1D<T> Ms = lTimesD.matrixTimesDiagonalVector(D_neg0pt5);
        
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
    
    
    DenseMatrix1D<T> scores;
    std::vector<int>* comp_mask_curr=comp_mask_values[0];
    int counter_eig_vector=0;
    int counter_comp_mask=0;
    struct IsoRank_Result ret_val;
    
    //run matching algorithms for each component
    for(int k=0;k<kron_prod.getNumberOfColumns();k++) {
        T* eigenvector=eigenValues[k];
        
        if(eigenvector!=NULL) {
            comp_mask_curr=comp_mask_values[k];
            scores= reshape(eigenvector,matrix_A.getNumberOfRows(),matrix_B.getNumberOfColumns(),*comp_mask_curr);
            DenseMatrix1D<T> scores_copy(scores);
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
                    DenseMatrix1D<T> perm_mat=getPermMatrix(assignment,matrix_A.getNumberOfRows(),matrix_A.getNumberOfRows());
                    DenseMatrix1D<T> product=perm_mat*matrix_A;
                    DenseMatrix1D<T> get_transpose=perm_mat.transpose();
                    DenseMatrix1D<T> final_mat=product*get_transpose;
                    DenseMatrix1D<T> ret_matrix= matrix_A-final_mat;
					
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
                    DenseMatrix1D<float> matrix_A2(b_size,b_size);
                    
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
                    
                    DenseMatrix1D<T> perm_mat=getPermMatrix(assignment,matrix_A.getNumberOfRows(),matrix_B.getNumberOfRows());
                    DenseMatrix1D<T> product=perm_mat*matrix_A2;	
                    DenseMatrix1D<T> get_transpose=perm_mat.transpose();
                    DenseMatrix1D<T> final_mat=product*get_transpose;
                    DenseMatrix1D<T> ret_matrix= matrix_A-final_mat;
                    
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

#endif
