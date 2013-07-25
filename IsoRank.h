//
//  IsoRank.h
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/14/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef Sparse_Matrix_IsoRank_h
#define Sparse_Matatorrix_IsoRank_h

#include "Matricies/Matrix2D.h"
#include "Matricies/SymSparseMatrix.h"
#include "Tarjan.h"
#include "util.h"
#include <vector>

const int GREEDY = 0;
const int CON_ENF_1 = 1;
const int CON_ENF_2 = 2;
const int CON_ENF_3 = 3;
const int CON_ENF_4 = 4;
const int NUM_OF_ISORANK_IT = 20;

struct IsoRank_Result
{
	int frob_norm;
    int assignment_length;
    int* assignments;
};

// void MPI_Send_IsoRank_Result (IsoRank_Result result, int dest, int tag)
// {
//     MPI_Send(&result.assignment_length, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
//     MPI_Send(result.assignments, result.assignment_length, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
//     MPI_Send(&result.frob_norm, 1, MPI_INT, dest, tag + 3, MPI_COMM_WORLD);
// }

// struct IsoRank_Result MPI_Recv_IsoRank_Result(int source, int tag, MPI_Status& stat)
// {
//     struct IsoRank_Result result;
//     MPI_Recv(&result.assignment_length, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
//     result.assignments = new int[result.assignment_length];
//     MPI_Recv(result.assignments ,result.assignment_length , MPI_INT, source, tag + 2, MPI_COMM_WORLD, &stat);
//     MPI_Recv(&result.frob_norm, 1, MPI_INT, source, tag + 3, MPI_COMM_WORLD, &stat);
//     return result;
// }

template <typename DT>
struct IsoRank_Result isoRank(SymSparseMatrix<DT>& matrix_A, SymSparseMatrix<DT>& matrix_B)
{

    // Degree distribution statistics
    
    SymSparseMatrix<DT> kron_prod = matrix_A.kron(matrix_B);
    
    //find componenets of the kronecker product and calculate the scores using eigenvector decomp
    vector<vertex*> * vertices = graph_con_com(kron_prod);
    double** eigenValues = new double*[kron_prod.getSize()];
    vector<vector<int>*> comp_mask_values (kron_prod.getSize());
    
    for(int i=0; i < kron_prod.getSize(); i++ )
    {
        vector<int>* comp_mask = component_mask(*vertices, i);
        comp_mask_values[i] = comp_mask;
        if (comp_mask == NULL)
        {
            eigenValues[i]=NULL;
            continue;
        }
        
        SparseMatrix<DT> L = kron_prod.getScatteredSelection(*comp_mask,*comp_mask);
        std::cout << L << std::endl;
//         
//         DT* sum = L->sum_rows();
//         DT* D_neg1  = new DT[L->getNumberOfRows()];
//         DT* D_0pt5  = new DT[L->getNumberOfRows()];
//         DT* D_neg0pt5  = new DT[L->getNumberOfRows()];
//         
//         for(int j=0; j < L->getNumberOfRows(); j++)
//         {
//             D_neg1[j] = 1.0/sum[j];
//             D_0pt5[j] = sqrt(sum[j]);
//             D_neg0pt5[j] = 1.0/D_0pt5[j];
//         }
//         
//         SparseMatrix<DT>* lTimesD = L->vec_times_mat(D_neg0pt5,L->getNumberOfRows());
//         SparseMatrix<DT>* Ms = lTimesD->mat_times_vec(D_neg0pt5, L->getNumberOfRows());
//         delete lTimesD;
//         if(!Ms->isSymmetric())
//         {
//             throw NotASymmetricMatrixException();
//         }
//         
//         //  cout<< "Ms Matrix" << endl << *Ms << endl;
//         
//         double* eigenVec=  Ms->getTopEigenVector();
//         double vecLength = 0;
//         for (int j=0; j < L->getNumberOfRows(); j++)
//         {
//             eigenVec[j] *= D_0pt5[j];
//             vecLength += pow(eigenVec[j],2);
//         }
//         
//         vecLength = sqrt(vecLength);
//         
//         int coef = 1;
//         if( eigenVec[0] < 0)
//         {
//             coef = -1;
//         }
//         
//         for (int j=0; j < L->getNumberOfRows(); j++)
//         {
//             eigenVec[j] = coef* (eigenVec[j]/vecLength);
//             //    cout << eigenVec[j] << ", ";
//         }
//         
//         eigenValues[i] = eigenVec;
//         
//         delete L;
//         delete [] sum;
//         delete [] D_neg1;
//         delete [] D_0pt5;
//         delete [] D_neg0pt5;
//         //delete M;
//         delete Ms;
//         
//         
//         //Algorithms
//         
    }
//     
//     
//     //reshape the eigenvector as a m by n matrix of floats and run
//     //node matching algorithms on it
//     SparseMatrix<double>* scores;
//     vector<int>* comp_mask_curr=comp_mask_values[0];
//     int counter_eig_vector=0;
//     int counter_comp_mask=0;
//     struct IsoRank_Result ret_val;
//     
//     for(int k=0;k<kron_prod->getNumberOfColumns();k++) {
//         double* eigenvector=eigenValues[k];
//         if(eigenvector!=NULL) {
//             comp_mask_curr=comp_mask_values[k];
//             scores= reshape(eigenvector,matrix_A.getNumberOfRows(),matrix_B.getNumberOfColumns(),*comp_mask_curr);
//             SparseMatrix<double>* scores_copy=new SparseMatrix<double>(*scores);
//             int * best_assignment=new int[matrix_A.getNumberOfRows()];
//             float best_frob_norm=DBL_MAX;
//             
//             for (int num_it = 0; num_it < NUM_OF_ISORANK_IT; num_it++){
//                 init_array(assignment,matrix_A.getNumberOfRows(),-1);
//                 *scores=*scores_copy;
//                 
//                 switch (matching_algorithm)
//                 {
//                     case 0:
//                         greedy_1(*scores,matrix_A,matrix_B,assignment);
//                         break;
//                     case 1:
//                         greedy_connectivity_1(*scores,matrix_A,matrix_B,assignment);
//                         break;
//                     case 2:
//                         greedy_connectivity_2(*scores,matrix_A,matrix_B,assignment);
//                         break;
//                     case 3:
//                         greedy_connectivity_3(*scores,matrix_A,matrix_B,assignment);
//                         break;
//                     case 4:
//                         greedy_connectivity_4(*scores,matrix_A,matrix_B,assignment);
//                         break;
//                     default:
//                         break;
//                 }
//                 
//                 //find the frobenius norm
//                 if(matrix_A.getNumberOfRows()>matrix_B.getNumberOfRows()){
//                     SparseMatrix<float> perm_mat=getPermMatrix(assignment,matrix_A.getNumberOfRows(),matrix_A.getNumberOfRows());
//                     SparseMatrix<float> product=perm_mat*matrix_A;
//                     SparseMatrix<float> get_transpose=perm_mat.transpose();
//                     SparseMatrix<float> final_mat=product*get_transpose;
//                     SparseMatrix<float> ret_matrix= matrix_A-final_mat;
//                     
//                     float frob_norm_hold=ret_matrix.getFrobNorm();
//                     
//                     if(frob_norm_hold<best_frob_norm){
//                         best_frob_norm=frob_norm_hold;
//                         best_assignment=assignment;
//                     }
//                 }
//                 else{
//                     int b_size=matrix_B.getNumberOfRows();
//                     int a_size=matrix_A.getNumberOfRows();
//                     SparseMatrix<float> matrix_A2(b_size,b_size);
//                     
//                     for(int r=0;r<b_size;r++){
//                         for(int s=0;s<b_size;s++){
//                             if(r<a_size&&s<a_size){
//                                 matrix_A2[r][s]=matrix_A[r][s];
//                             }
//                             else if(r==s){
//                                 matrix_A2[r][s]=1;
//                             }
//                             else{
//                                 matrix_A2[r][s]=0;
//                             }
//                             
//                         }
//                     }
//                     
//                     SparseMatrix<float> perm_mat=getPermMatrix(assignment,matrix_A.getNumberOfRows(),matrix_B.getNumberOfRows());
//                     SparseMatrix<float> product=perm_mat*matrix_A2;	
//                     SparseMatrix<float> get_transpose=perm_mat.transpose();
//                     SparseMatrix<float> final_mat=product*get_transpose;
//                     SparseMatrix<float> ret_matrix= matrix_A-final_mat;
//                     
//                     float frob_norm_hold=ret_matrix.getFrobNorm(); 
//                     
//                     if(frob_norm_hold<best_frob_norm){
//                         best_frob_norm=frob_norm_hold;
//                         best_assignment=assignment;
//                     }
//                 }
//             }
//             
//             ret_val.frob_norm=best_frob_norm;
//             ret_val.assignments=best_assignment;
//             ret_val.assignment_length=matrix_A.getNumberOfRows();
//             delete []best_assignment;
//             
//         }
//         return ret_val;
//     }
//     
//     //delete memory
//     
//     for(int k=0;k<kron_prod->getNumberOfColumns();k++)
//     {
//         delete [] eigenValues[k];
//     }
//     
//     delete scores;
//     typename vector<vector<int>*>::iterator i;
//     
//     
//     for ( i = comp_mask_values.begin() ; i < comp_mask_values.end(); ++i ){
//         delete  * i;
//     }
//     
//     
//     typename vector<vertex*>::iterator it;
//     for ( it = vertices->begin() ; it < vertices->end(); ++it )
//     {
//         delete * it;
//     }
//     
//     delete vertices;
//     delete kron_prod;
//     delete []eigenValues;
    
}

#endif
