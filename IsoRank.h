//
//  IsoRank.h
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/14/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef Sparse_Matrix_IsoRank_h
#define Sparse_Matrix_IsoRank_h

#include "SparseMatrix.h"
#include "Tarjan.h"
#include "util.h"
#include  "greedy_algorithms.h"
#include <vector>

int GREEDY = 0;
int CON_ENF_1 = 1;
int CON_ENF_2 = 2;
int CON_ENF_3 = 3;
int CON_ENF_4 = 4;



template <typename DT>
void isoRank(SparseMatrix<DT>& matrix_A, SparseMatrix<DT>& matrix_B, int matching_algorithm)
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
    
    
     SparseMatrix<DT>* kron_prod = matrix_A.kron(matrix_B);
    
   
    vector<vertex*> * vertices = graph_con_com(kron_prod, kron_prod->getNumberOfColumns());
    double* eigenValues[kron_prod->getNumberOfColumns()];
    vector<vector<int>*> comp_mask_values (kron_prod->getNumberOfColumns()); 

    for(int i=0; i < kron_prod->getNumberOfColumns(); i++ )
    {
        vector<int>* comp_mask = component_mask(*vertices, i);
        comp_mask_values[i] = comp_mask;
        if (comp_mask == NULL)
        {
	  eigenValues[i]=NULL;
            continue;
        }
    
        SparseMatrix<DT>* L = kron_prod->getScatteredSelection(*comp_mask,*comp_mask);
        
        DT* sum = L->sum_rows();
        DT* D_neg1  = new DT[L->getNumberOfRows()];
        DT* D_0pt5  = new DT[L->getNumberOfRows()];
        DT* D_neg0pt5  = new DT[L->getNumberOfRows()];
        
        for(int j=0; j < L->getNumberOfRows(); j++)
        {
            D_neg1[j] = 1.0/sum[j];
            D_0pt5[j] = sqrt(sum[j]);
            D_neg0pt5[j] = 1.0/D_0pt5[j];
        }
    
        //SparseMatrix<DT>* M = L->vec_times_mat(D_neg1,L->getNumberOfRows() );
		SparseMatrix<DT>* lTimesD = L->vec_times_mat(D_neg0pt5,L->getNumberOfRows());
        SparseMatrix<DT>* Ms = lTimesD->mat_times_vec(D_neg0pt5, L->getNumberOfRows());
		delete lTimesD;
        if(!Ms->isSymmetric())
        {
            throw NotASymmetricMatrixException();
        }
    
	//  cout<< "Ms Matrix" << endl << *Ms << endl;

        double* eigenVec=  Ms->getTopEigenVector();
        double vecLength = 0;
        for (int j=0; j < L->getNumberOfRows(); j++)
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
        
        for (int j=0; j < L->getNumberOfRows(); j++)
        {
            eigenVec[j] = coef* (eigenVec[j]/vecLength);
	    //    cout << eigenVec[j] << ", ";
        }
        
        eigenValues[i] = eigenVec;
        
        delete L;
        delete [] sum;
        delete [] D_neg1;
        delete [] D_0pt5;
        delete [] D_neg0pt5;
        //delete M;
        delete Ms;
    

    //Algorithms

    }


    SparseMatrix<DT>* scores= new SparseMatrix<DT>(matrix_A.getNumberOfRows(),matrix_B.getNumberOfColumns());
    vector<int>* comp_mask_curr=comp_mask_values[0];
    int counter_eig_vector=0;
    int counter_comp_mask=0;
    
    
    for(int k=0;k<kron_prod->getNumberOfColumns();k++){
        double* eigenvector=eigenValues[k];
        if(eigenvector!=NULL) {
            for(int j=0;j<scores->getNumberOfRows();j++){
                for(int i=0;i<scores->getNumberOfColumns();i++){
                    if((*comp_mask_curr)[counter_comp_mask]==1){
                        (*scores)[j][i]=eigenvector[counter_eig_vector];
                        counter_eig_vector++;
                    }
                    else{
                        (*scores)[j][i]=0;
                    }
                    counter_comp_mask++;
                }
            
            
            
            int* assignment = new int[matrix_A.getNumberOfRows()];
            
            switch (matching_algorithm)
            {
                case 0:
                    greedy_1(*scores,assignment);
                    break;
                case 1:
                     greedy_connectivity_1(*scores,matrix_A,matrix_B,assignment);
                    break;
                case 2:
                     greedy_connectivity_2(*scores,matrix_A,matrix_B,assignment);
                    break;
                case 3:
                     greedy_connectivity_3(*scores,matrix_A,matrix_B,assignment);
                    break;
                case 4:
                     //greedy_connectivity_4(*scores,matrix_A,matrix_B,assignment);
                    break;
                default:
                    break;
            }
           
            
            // int* assignment= greedy_1(*scores);
            for(int i=0;i<matrix_A.getNumberOfRows();i++){
                printf(" graph1: %d graph2: %d\n",i,assignment[i]);
            }
			delete [] assignment;
        }
}
    }
    

  //comp_mask_values[i]
    //eigenValues[i]
	 for(int k=0;k<kron_prod->getNumberOfColumns();k++)
	{
        delete [] eigenValues[k];
	}
	delete scores;
	 typename vector<vector<int>*>::iterator i;
    for ( i = comp_mask_values.begin() ; i < comp_mask_values.end(); ++i )
    {
        delete  * i;
    }


    typename vector<vertex*>::iterator it;
    for ( it = vertices->begin() ; it < vertices->end(); ++it )
    {
        delete * it;
    }
    delete vertices;
    delete kron_prod;
	
    
}

#endif
