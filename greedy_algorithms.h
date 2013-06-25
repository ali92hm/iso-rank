#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "SparseMatrix.h"


#ifndef greedy_algorithm_h
#define greedy_algorithm_h

template <typename DT>
SparseMatrix<DT>* reshape(double* &eigenvector,const int rows,const int cols, vector<int>* &comp_mask){
  SparseMatrix<DT>* matrix= new SparseMatrix<DT>(rows,cols);
  int counter_eig_vector=0;
  int counter_comp_mask=0;
  for(int j=0;j<matrix->getNumberOfColumns();j++){
    for(int i=0;i<matrix->getNumberOfRows();i++){
      if((*comp_mask)[counter_comp_mask]==1){
	(*matrix)[i][j]=eigenvector[counter_eig_vector];
	counter_eig_vector++;
      }
      else{
	(*matrix)[i][j]=0;
      }
      counter_comp_mask++;
    }
  }

  return matrix;
}




template <typename DT>
int* greedy_1(SparseMatrix<DT>& matches){
  DT total_score=0;
  int graph1_nodes=matches.getNumberOfRows();
  int graph2_nodes=matches.getNumberOfColumns();

  int num_of_nodes=min(graph1_nodes,graph2_nodes);
  int max_value;
  int* assignment= (int*)malloc(sizeof(int)*graph1_nodes);
  int row,col;
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
  }

  for(int i=0;i<min(graph1_nodes,graph2_nodes);i++){
    return_max(matches,&total_score,&row,&col);
    invalidate(row,col,matches);
    assignment[row]=col;
  }
  
  return assignment;
}


template<typename DT>
int* greedy_connectivity_1(SparseMatrix<DT>& matches, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2,int* assignment){
  
  DT total_score=0;
  int graph1_nodes=matches.getNumberOfRows();
  int graph2_nodes=matches.getNumberOfColumns();

  int max_value;
  
  int row,col;
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
  }
 
  for(int i=0;i<min(graph1_nodes,graph2_nodes);i++){
    return_max(matches,&total_score,&row,&col);
    assignment[row]=col;
    invalidate(row,col,matches);
    neighbor_enforcement(&row,&col, graph1,graph2,matches);
    total_score+=max_value;
  }

  return assignment;
}






template<typename DT>
void neighbor_enforcement(int* row_index,int* col_index, SparseMatrix<float>& graph1,SparseMatrix<float>& graph2, SparseMatrix<DT>& matches){
    for(int i=0;i<graph1.getNumberOfColumns();i++){
      if(graph1[*row_index][i]==1){
	for(int j=0;j<graph2.getNumberOfRows();j++){
	  if(graph2[j][*col_index]==0){
	    matches[i][j]=-DBL_MAX;
	  }
	}
      }
    }
  
    for(int i=0;i<graph2.getNumberOfColumns();i++){
      if(graph2[*col_index][i]==1){
	for(int j=0;j<graph1.getNumberOfRows();j++){
	  if(graph1[j][*row_index]==0){
	    matches[j][i]==-DBL_MAX;
	  }
	}
      }
    }
}



template<typename DT>
int* greedy_connectivity_2(SparseMatrix<DT>& matches, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2,int* assignment){
  DT max_tol=pow(10,-6);
  DT score=0,prev_score=0;
  
  int graph1_nodes=graph1.getNumberOfColumns();
  int graph2_nodes=graph2.getNumberOfColumns();

  int assignment2[graph1_nodes];
  int row,col,max;
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
    assignment2[i]=0;
  }

  SparseMatrix<DT>* active_matches=new SparseMatrix<DT>(matches);
  DT* idxarr;
  vector<int> assigned_G1;
  int size=0,random_id,vector_size=0,curr_row;
 
srand(time(NULL));
  while(sum_array(assignment2,graph1_nodes)<min(graph1_nodes,graph2_nodes) && (max = return_max(matches,&prev_score,&row,&col))>-DBL_MAX){
      
     
      idxarr=find_values(matches, prev_score-max_tol,&size);
      prev_score=0;

      random_id=rand()%size;
      printf("random id: %d size: %d\n",random_id,size);
      max= get_Max(&matches,random_id,max-max_tol,&row,&col);
      score+=max;
      assignment[row]=col;
      assignment2[row]=1;
      
      assigned_G1.push_back(row);
      vector_size++;
      invalidate(row, col, matches);
      
      set_to_min(*active_matches);
      int curr_col;

      for(int i=0;i<vector_size;i++) {
	curr_row=assigned_G1[i];
	curr_col=assignment[curr_row];
	vector<int>* neigh_1= graph1.getNeighbors(curr_row);
	vector<int>* neigh_2= graph2.getNeighbors(curr_col);
	set_matrix_values(*active_matches,matches,*neigh_1,*neigh_2);
      }
  }

}




template<typename DT>
int* greedy_connectivity_3(SparseMatrix<DT>& matches, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2,int* assignment){
  DT total_score=0;
  int row,col;
  int graph1_nodes=graph1.getNumberOfRows();
  int graph2_nodes=graph2.getNumberOfRows();
  assignment= malloc(sizeof(int)* graph1.getNumberOfRows());
  int assignment2[graph1.getNumberOfRows()];
  int assignment_G1[graph1.getNumberOfRows()];
  int assignment_G2[graph1.getNumberOfRows()];;
  
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
    assignment_G1[i]=0;
    assignment_G2[i]=0;
  }
  vector<int>* neigh_1;
  vector<int>* neigh_2;

  while(sum_array(assignment2,graph1_nodes)<min(graph1_nodes,graph2_nodes)){
    return_max(matches, &total_score,&row,&col);
    assignment[row]=col;
    assignment2[row]=1;
    invalidate(row,col,matches);
    neigh_1= graph1.getNeighbors(row);
    neigh_2= graph2.getNeighbors(col);
    
    for(int i=0;i<neigh_1->size();i++){
      if(assignment_G1[neigh_1[i]]!=1){
	neigh_1->erase(i);
      }

    }

    for(int i=0;i<neigh_2->size();i++){
       if(assignment_G2[neigh_2[i]]==1){
	 neigh_2->erase(i);
      }

    }
    


  }


  
  

}


template <typename DT>
void set_to_min(SparseMatrix<DT>& matrix){
 
  for(int i=0;i<matrix.getNumberOfRows();i++){
    for(int j=0;j<matrix.getNumberOfColumns();j++){
      matrix[i][j]=-DBL_MAX;
    }
  }

}

template<typename DT>
void set_matrix_values(SparseMatrix<DT>& matrix1,SparseMatrix<DT>& matrix2,  vector<int>& rows, vector<int>& cols){
  int r;
  int c;

  for(int i=0;i<rows.size();i++){
     r=rows[i];
    for(int j=0;j<cols.size();j++){
      c=cols[j];
      matrix1[r][c]=matrix2[r][c];
    }
  }


}


template <typename DT>
DT* find_values(SparseMatrix<DT>& matches2,DT value,int *size){
  *size=0;

  for(int i=0; i<matches2.getNumberOfRows();i++){
    for(int j=0;j<matches2.getNumberOfColumns();j++){
      if(matches2[i][j]>=value){
	*size=*size+1;
      }
    }
  }

  DT* retarr;
  
  retarr=(DT*)malloc(sizeof(DT)*(*size));
  int counter=0;
  
  
  for(int i=0;i<matches2.getNumberOfRows();i++){
    for(int j=0;j<matches2.getNumberOfColumns();j++){
      if(matches2[i][j]>=value){
	 retarr[counter]=matches2[i][j];
	  counter++;
	}
    }
    }
 
  return retarr;
} 

template <typename DT>
DT get_Max(SparseMatrix<DT> *matches, int id, DT val,int *row, int *col){
  int count=0;

  for(int i=0; i<matches->getNumberOfRows();i++){
    for(int j=0;j<matches->getNumberOfColumns();j++){
      if((*matches)[i][j]>=val){
	count++;
	if(count==id){
	  *row=i;
	  *col=j;
	  return (*matches)[i][j];
	}
      }
    }
  }
}




template <typename DT>
int return_max(SparseMatrix<DT>& matches, DT* total_score,int* max_row,int* max_col){
  DT max_so_far=matches[0][0];
  int max_so_far_count=1;
  *max_row=0; 
  *max_col=0;

  for(int i=0; i<matches.getNumberOfRows();i++){
    for(int j=0;j<matches.getNumberOfColumns();j++){
      if(matches[i][j]==max_so_far){
	max_so_far_count++;
      }
      else if(matches[i][j]>max_so_far){
	max_so_far=matches[i][j];
	max_so_far_count=1;
      }
      
    }
  }
  srand (time(NULL));

  int counter=0;
  int random_number= rand()%max_so_far_count+1;
  int set=0;

  printf("random_number: %d and max so far count:%d\n",random_number, max_so_far_count);

  for(int i=0; i<matches.getNumberOfRows();i++){
    for(int j=0;j<matches.getNumberOfColumns();j++){
       if(matches[i][j]==max_so_far){
	 counter++;
	 if(counter==random_number){
	   *max_row=i;
	   *max_col=j;
	   set=1;
	   break;
	 }
      }

     
    }
    if(set==1)
      break;

  }

 

   *total_score+=max_so_far;

  return *max_col;
}

template <typename DT>
void invalidate(int row, int col, SparseMatrix<DT>& matches){
  for(int j=0;j<matches.getNumberOfColumns();j++){
    matches[row][j]=-DBL_MAX;
  }
  for(int j=0;j<matches.getNumberOfRows();j++){
    matches[j][col]=-DBL_MAX;
  }

}
#endif
