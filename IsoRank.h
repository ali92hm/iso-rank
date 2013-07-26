//
//  IsoRank.h
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/14/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef Sparse_Matrix_IsoRank_h
#define Sparse_Matatorrix_IsoRank_h

#include "Matricies/DenseMatrix.h"
#include "Tarjan.h"
#include "util.h"
#include "greedy_algorithms.h"
#include <vector>

int GREEDY = 0;
int CON_ENF_1 = 1;
int CON_ENF_2 = 2;
int CON_ENF_3 = 3;
int CON_ENF_4 = 4;

const int NUM_OF_ISORANK_IT = 20;

struct IsoRank_Result
{
    
	int frob_norm;
    int assignment_length;
	int* assignments;
    
    
};

void MPI_Send_IsoRank_Result (IsoRank_Result result, int dest, int tag)
{
	MPI_Send(&result.assignment_length, 1, MPI_INT, dest, tag + 1, MPI_COMM_WORLD);
	MPI_Send(result.assignments, result.assignment_length, MPI_INT, dest, tag + 2, MPI_COMM_WORLD);
	MPI_Send(&result.frob_norm, 1, MPI_INT, dest, tag + 3, MPI_COMM_WORLD);
}

struct IsoRank_Result MPI_Recv_IsoRank_Result(int source, int tag, MPI_Status& stat)
{
// 	std::cout << "CALL to MPI_Recv_IsoRank_Result "<< std::endl;
	struct IsoRank_Result result;
	MPI_Recv(&result.assignment_length, 1, MPI_INT, source, tag + 1, MPI_COMM_WORLD, &stat);
	result.assignments = new int[result.assignment_length];
	MPI_Recv(result.assignments ,result.assignment_length , MPI_INT, source, tag + 2, MPI_COMM_WORLD, &stat);
	MPI_Recv(&result.frob_norm, 1, MPI_INT, source, tag + 3, MPI_COMM_WORLD, &stat);
	return result;
}

template <typename DT>
struct IsoRank_Result isoRank(DenseMatrix<DT>& matrix_A, DenseMatrix<DT>& matrix_B, int matching_algorithm)
{
    if (!matrix_A.isSquare() || !matrix_B.isSquare())
    {
        throw NotASquareMatrixException();
    }
    
    if (!matrix_A.isSymmetric() || !matrix_B.isSymmetric())
    {
        throw NotASymmetricMatrixException();
    }
    
    // Degree distribution statistics
    
    
    DenseMatrix<DT> kron_prod = matrix_A.kron(matrix_B);
    std::vector<vertex*> vertices = graph_con_com(kron_prod);
    DT** eigenValues = new DT*[kron_prod.getNumberOfColumns()];
    vector<vector<int>*> comp_mask_values (kron_prod.getNumberOfColumns());
    
    for(int i=0; i < kron_prod.getNumberOfColumns(); i++ )
    {
        vector<int>* comp_mask = component_mask(vertices, i);
        comp_mask_values[i] = comp_mask;
        if (comp_mask == NULL)
        {
	  		eigenValues[i]=NULL;
            continue;
        }
        
        DenseMatrix<DT> L = kron_prod.getScatteredSelection(*comp_mask,*comp_mask);
        
        std::vector<DT> sum = L.getSumOfRows();
        std::vector<DT> D_neg1(sum.size());
        std::vector<DT> D_0pt5(sum.size());
        std::vector<DT> D_neg0pt5(sum.size());
        
        for(int j=0; j < sum.size(); j++)
        {
            D_neg1[j] = 1.0/sum[j];
            D_0pt5[j] = sqrt(sum[j]);
            D_neg0pt5[j] = 1.0/D_0pt5[j];
        }
        
		DenseMatrix<DT> lTimesD = L.diagonalVectorTimesMatrix(D_neg0pt5);
        DenseMatrix<DT> Ms = lTimesD.matrixTimesDiagonalVector(D_neg0pt5);
        
        if(!Ms.isSymmetric())
        {
            throw NotASymmetricMatrixException();
        }
        
        DT* eigenVec=  Ms.getTopEigenVector();
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
    
    
    DenseMatrix<DT> scores;
    vector<int>* comp_mask_curr=comp_mask_values[0];
    int counter_eig_vector=0;
    int counter_comp_mask=0;
    struct IsoRank_Result ret_val;
    //int * best_assignment = new int[matrix_A.getNumberOfRows()];
  
    for(int k=0;k<kron_prod.getNumberOfColumns();k++) {
        DT* eigenvector=eigenValues[k];
        if(eigenvector!=NULL) {
      comp_mask_curr=comp_mask_values[k];
      scores= reshape(eigenvector,matrix_A.getNumberOfRows(),matrix_B.getNumberOfColumns(),*comp_mask_curr);
      DenseMatrix<DT> scores_copy(scores);
      int * best_assignment/*=new int[matrix_A.getNumberOfRows()]*/;
      float best_frob_norm=DBL_MAX;
            
      for (int num_it = 0; num_it < NUM_OF_ISORANK_IT; num_it++){
        //init_array(assignment,matrix_A.getNumberOfRows(),-1);
        scores=scores_copy;
        int* assignment= new int[matrix_A.getNumberOfRows()];
		init_array(assignment,matrix_A.getNumberOfRows(),-1);
        switch (matching_algorithm)
          {
          case 0:
        greedy_1(scores,matrix_A,matrix_B,assignment);
        break;
          case 1:
        greedy_connectivity_1(scores,matrix_A,matrix_B,assignment);
        break;
          case 2:
        greedy_connectivity_2(scores,matrix_A,matrix_B,assignment);
        break;
          case 3:
        greedy_connectivity_3(scores,matrix_A,matrix_B,assignment);
        break;
          case 4:
        greedy_connectivity_4(scores,matrix_A,matrix_B,assignment);
        break;
          default:
        break;
          }
                
        //find the frobenius norm               
       if(matrix_A.getNumberOfRows()>matrix_B.getNumberOfRows()){
	      DenseMatrix<float> perm_mat=getPermMatrix(assignment,matrix_A.getNumberOfRows(),matrix_A.getNumberOfRows());
	      DenseMatrix<float> product=perm_mat*matrix_A;	
	      DenseMatrix<float> get_transpose=perm_mat.transpose();	
	      DenseMatrix<float> final_mat=product*get_transpose;
	      DenseMatrix<float> ret_matrix= matrix_A-final_mat;
					
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
	    else{
	      int b_size=matrix_B.getNumberOfRows();
	      int a_size=matrix_A.getNumberOfRows();
	      DenseMatrix<float> matrix_A2(b_size,b_size);
				  
	      for(int r=0;r<b_size;r++){
			for(int s=0;s<b_size;s++){
			if(r<a_size&&s<a_size){
		 	   matrix_A2[r][s]=matrix_A[r][s];
		 	 }
		 	 else if(r==s){
		 	   matrix_A2[r][s]=1;
			  }
		 	 else{
		 	   matrix_A2[r][s]=0;
		 	 }	      
			}
	      }

	      DenseMatrix<float> perm_mat=getPermMatrix(assignment,matrix_A.getNumberOfRows(),matrix_B.getNumberOfRows());
	      DenseMatrix<float> product=perm_mat*matrix_A2;	
	      DenseMatrix<float> get_transpose=perm_mat.transpose();
	      DenseMatrix<float> final_mat=product*get_transpose;
	      DenseMatrix<float> ret_matrix= matrix_A-final_mat;

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