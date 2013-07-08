#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include "SparseMatrix.h"

#ifndef greedy_algorithm_h
#define greedy_algorithm_h



struct coordinate_pair{
  int row;
  int col;
};

vector<int>* intersect(int*, int, struct coordinate_pair**,int);
void match_rest(int*, SparseMatrix<float>&, SparseMatrix<float>&);
int* get_valid_entries(SparseMatrix<float>, int*,int,int*);



/*
 * performs greedy algorithm on matrix of nodal pairs matched and 
 * returns a matching between nodes of graph1 and graph2
 * @pram: matrix indicating the scores of nodal pairings
 */
template <typename DT>
int* greedy_1(SparseMatrix<DT>& matches,int* assignment){
  DT total_score=0;
  int graph1_nodes=matches.getNumberOfRows();
  int graph2_nodes=matches.getNumberOfColumns();

  int num_of_nodes=min(graph1_nodes,graph2_nodes);
  int max_value;
  int row,col;
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
  }

  for(int i=0;i<min(graph1_nodes,graph2_nodes);i++){
    return_max(matches,&total_score,&row,&col);
    invalidate(row,col,matches);
    assignment[row]=col;
  }
  printf("score of matching: %f\n",total_score);
  return assignment;
}


/*
 * peforms a greedy algorithm to choose the best nodal pairs for matching 
 * enforces connectivity: if i<->j then neigh(i)<->neigh(j) where <-> indicates a matching
 * @pram: matrix indicating scores for nodal pairs 
 * @pram: adjacency matrix for graph1
 * @pram: adjacency matrix for graph2
 * @pram: pointer to the array that indicates the best matching
 */
template<typename DT>
int* greedy_connectivity_1(SparseMatrix<DT>& matches, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2,int* assignment){
  
  DT total_score=0;
  int graph1_nodes=matches.getNumberOfRows();
  int graph2_nodes=matches.getNumberOfColumns();

  DT max_value;
  
  int row,col;
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
  }
 
  for(int i=0;i<min(graph1_nodes,graph2_nodes);i++){
    return_max(matches,&total_score,&row,&col);
    assignment[row]=col;
    invalidate(row,col,matches);
    neighbor_enforcement(&row,&col, graph1,graph2,matches);
  }

  printf("score of matching: %f\n",total_score);
  return assignment;
}



/*
 * performs a greedy matching and enforces connectivity by proceeding outwards radially
 * @pram:matrix indicating scores for nodal pairs 
 * @pram: adjacency matrix for graph1
 * @pram: adjacency matrix for graph2
 * @pram: pointer to the array that indicates the best matching 
 */
template<typename DT>
void greedy_connectivity_2(SparseMatrix<DT>& matches, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2,int* assignment){
  DT max_tol=pow(10,-6),max;
  DT score=0,prev_score=0,final_score=0;
  
  int graph1_nodes=graph1.getNumberOfColumns();
  int graph2_nodes=graph2.getNumberOfColumns();
  int assignment2[graph1_nodes];
  int row,col;
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
    assignment2[i]=0;
  }

  SparseMatrix<DT>* active_matches=new SparseMatrix<DT>(matches);
  DT* idxarr;
  vector<int> assigned_G1;
  int size=0,random_id,vector_size=0,curr_row;
 
srand(time(NULL));
 while(sum_array(assignment2,graph1_nodes)<min(graph1_nodes,graph2_nodes)&& return_max(*active_matches,&score,&row,&col)>-1)
{

    if(all_inf(*active_matches)){
        match_rest(assignment,graph1,graph2); 
        printf("score of matching: %f\n", final_score);
	return;
      }
 
    // idxarr=find_values(*active_matches, score-prev_score-max_tol,&size);
      random_id=rand()%size+1;
      max= get_Max(active_matches,random_id,score-prev_score-max_tol,&row,&col);
    
      final_score+=max;
      assignment[row]=col;
      assignment2[row]=1;
      
      assigned_G1.push_back(row);
      vector_size++;
      invalidate(row, col, matches);
      set_to_min(*active_matches);
      int curr_col;

      for(int i=0;i<vector_size;i++) 
      {
	curr_row=assigned_G1[i];
	curr_col=assignment[curr_row];
	vector<int>* neigh_1= graph1.getNeighbors(curr_row);
	vector<int>* neigh_2= graph2.getNeighbors(curr_col);
	set_matrix_values(*active_matches,matches,*neigh_1,*neigh_2);
      }
      prev_score=score;
  }

 printf("score of matching: %f\n", final_score);
}



/*
 * performs a greedy matching and enforces connectivity by proceeding outwards radially
 * @pram:matrix indicating scores for nodal pairs 
 * @pram: adjacency matrix for graph1
 * @pram: adjacency matrix for graph2
 * @pram: pointer to the array that indicates the best matching 
 */
template<typename DT>
void greedy_connectivity_3(SparseMatrix<DT>& matches, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2,int* assignment){
  DT total_score=0;
  DT final_score=0;
  int row,col;
  int graph1_nodes=graph1.getNumberOfRows();
  int graph2_nodes=graph2.getNumberOfRows();
  int assignment2[graph1.getNumberOfRows()];
  int assignment_G1[graph1.getNumberOfRows()];
  int assignment_G2[graph2.getNumberOfRows()];;
  SparseMatrix<DT>* local_matches=new SparseMatrix<DT>(matches);
  

  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
    assignment2[i]=0;
    assignment_G1[i]=0;
    assignment_G2[i]=0;
  }

  vector<int>* neigh_1;
  vector<int>* neigh_2;

  while(sum_array(assignment2,graph1_nodes)<min(graph1_nodes,graph2_nodes)){
    return_max(matches, &final_score,&row,&col);
    assignment[row]=col;
    assignment2[row]=1;
    assignment_G1[row]=1;
    assignment_G2[col]=1;
    invalidate(row,col,matches);
    neigh_1= graph1.getNeighbors(row);
    neigh_2= graph2.getNeighbors(col);
    int hold;
     for(int i=0;i<neigh_1->size();i++){
      hold=(*neigh_1)[i];
      if(assignment_G1[hold]==1){
	neigh_1->erase(neigh_1->begin()+i);
      }
    }

    for(int j=0;j<neigh_2->size();j++){
      hold=(*neigh_2)[j];
      if(assignment_G2[hold]==1){
	neigh_2->erase(neigh_2->begin()+j);
      }
    }

    set_to_min(*local_matches);
    set_matrix_values(*local_matches,matches,*neigh_1,*neigh_2);
    int row_inside= row;
    int col_inside=col;
     if(all_inf(matches)){
       match_rest(assignment,graph1,graph2);
       return;
     }

    
    for(int i=0;i<min(neigh_1->size(),neigh_2->size())&&!all_inf(*local_matches);i++){
      return_max(*local_matches,&final_score,&row,&col);
      assignment[row]=col;
      assignment2[row]=1;
      assignment_G1[row]=1;
      assignment_G2[col]=1;
      invalidate(row,col,*local_matches);
      invalidate(row,col,matches);

        if(all_inf(matches)){
	  match_rest(assignment,graph1,graph2);
	  return;
	}

    }
  }
  return;
}




template <typename DT>
void greedy_connectivity_4(SparseMatrix<DT>& matches, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2,int* assignment,int num_iterations){
  DT optScore=0;
  vector<DT> broken1;
  vector<DT> broken2;
  DT idx_b=0;
  DT final_score=0;
  int itr=0;
  int* add_order=new int[graph1.getNumberOfRows()];
  int* ass=new int[graph1.getNumberOfRows()];
  int add_idx=1;
  DT score=0;
  DT max_tol=pow(10,-6);
  int row,col,size=0;
  int assigned_G1[graph1.getNumberOfRows()];
  int assigned_G2[graph2.getNumberOfRows()];
  vector<int>* neigh_1;
  vector<int>* neigh_2;

 for(int j=0;j<graph1.getNumberOfRows();j++){
    add_order[j]=0;
    ass[j]=0;
    assigned_G1[j]=-1;
    assigned_G2[j]=-1;
  }

  
 return_max(matches, &score,&row,&col);
 DT* idx_array=find_values(matches,score - max_tol,&size);
 int random_id=rand()%size+1; 
 DT max= get_Max(&matches,random_id,score-max_tol,&row,&col);
 final_score+=max;
 assignment[row]=col;
 assigned_G1[row]=1;
 assigned_G2[col]=1;
 invalidate(row,col,matches);
 add_order[0]=row;
 neigh_1= graph1.getNeighbors(row);
 neigh_2= graph2.getNeighbors(col);

 int hold;
     for(int i=0;i<neigh_1->size();i++){
      hold=(*neigh_1)[i];
      if(assigned_G1[hold]==1){
	neigh_1->erase(neigh_1->begin()+i);
      }
    }

    for(int j=0;j<neigh_2->size();j++){
      hold=(*neigh_2)[j];
      if(assigned_G2[hold]==1){
	neigh_2->erase(neigh_2->begin()+j);
      }
    }

   SparseMatrix<DT>* matches_local=new SparseMatrix<DT>(matches.getNumberOfRows(),matches.getNumberOfColumns());
    set_matrix_values(*matches_local, matches, *neigh_1, *neigh_2);
    score=0;
    size=0;
    return_max(*matches_local, &score,&row,&col);
    idx_array=find_values(*matches_local,score-max_tol,&size);
    random_id=rand()%size+1;
    max=get_Max(matches_local,random_id,score-max_tol,&row,&col);
    final_score+=max;
    assignment[row]=col;
    assigned_G1[row]=1;
    assigned_G2[col]=1;
    invalidate(row,col,*matches_local);
    invalidate(row,col,matches);
    add_order[1]=row;
    
   for(int i=0;i<neigh_1->size();i++){
      hold=(*neigh_1)[i];
      if(assigned_G1[hold]==1){
	neigh_1->erase(neigh_1->begin()+i);
      }
    }

    for(int j=0;j<neigh_2->size();j++){
      hold=(*neigh_2)[j];
      if(assigned_G2[hold]==1){
	neigh_2->erase(neigh_2->begin()+j);
      }
    }
    int best_row,best_col;
    
    while(add_order[graph1.getNumberOfRows()-1]==0){
      for(int i=0; i<min(neigh_1->size(),neigh_2->size());i++) {
	score=0;
	return_max(*matches_local,&score,&row,&col);
	idx_array=find_values(*matches_local,score-max_tol,&size);
	int *rows_cols_size=0;
	int valid_entries_size=0,valid_entries2_size=0;
	struct coordinate_pair **rows_cols=find_all_values(*matches_local,idx_array,size,rows_cols_size);
	if(size>1) {
	  int* valid_entries= get_valid_entries(graph1,assignment,size,&valid_entries_size);
	  vector<int>* prev_assigned = intersect(valid_entries,valid_entries_size,rows_cols,size);
	  vector<int>* g1c_count;
	  int g1c_count_counter=0;
	  for(int k=0;k<prev_assigned->size();k++){
	    int id=(*prev_assigned)[k];
	    int sum=0;

	    for(int j=0;j<valid_entries_size;j++){
	      if(valid_entries[j]==id){
		sum++;
	      }
	    }

	    g1c_count->insert(g1c_count->begin()+g1c_count_counter,sum);
	    g1c_count_counter++;
	  }
	  
	  vector<int> *max_g1c=vector_max(g1c_count);
	  int rand_number = rand()%(max_g1c->size());
	  best_row=(*max_g1c)[rand_number];
	  vector<int>* best_cols= new vector<int>();
	  int best_cols_counter=0;
	  for(int i=0;i<*rows_cols_size;i++){
	    struct coordinate_pair *c_pair=rows_cols[i];
	    if(c_pair->row==best_row){
	      best_cols->insert(best_cols->begin()+best_cols_counter,c_pair->col);
	      best_cols_counter++;
	    }
	  }

	  int* valid_entries2= get_valid_entries(graph2,assignment,size,&valid_entries2_size);

	  vector<int>* g2c_count;
	  int g2c_count_counter=0;
	  for(int k=0;k<valid_entries2_size;k++){
	    int id=(*prev_assigned)[k];
	    int sum=0;
	//  vector<int> prev_assigned = intersect(valid_entries,valid_entries_size,rows_cols,size);

	    for(int j=0;j<valid_entries2_size;j++){
	      if(valid_entries2[j]==id){
		sum++;
	      }
	    }

	    g2c_count->insert(g2c_count->begin()+g2c_count_counter,sum);
	    g2c_count_counter++;
	  }

	  vector<int> *max_g2c=vector_max(g2c_count);
	   rand_number = rand()%(max_g2c->size());
	   best_col=(*max_g2c)[rand_number];
	}
	else{
	  struct coordinate_pair *cp=rows_cols[0];
	  best_row=cp->row;
	  best_col=cp->col;
	}
	assignment[best_row]=best_col;
	assigned_G1[best_row]=1;
	assigned_G2[best_col]=1;
	DT max_matches=matches[best_row][best_col];
	invalidate(best_row,best_col,*matches_local);
	invalidate(best_row,best_col,matches);
	final_score+=max_matches;
	
      }


    }
}







/*
 * returns a SparseMatrix Object which is a reshaped eigenvector
 * @pram: a pointer to an array of doubles which represents the eigenvector
 * @pram: number of rows in the matrix returned
 * @pram: number of columns in the matrix returned
 * @pram: component mask vector indicating which nodes are present in current component
 */
template <typename DT>
SparseMatrix<DT>* reshape(DT* eigenvector,const int rows,const int cols, vector<int> &comp_mask){

  SparseMatrix<DT>* matrix= new SparseMatrix<DT>(rows,cols);
  int counter_eig_vector=0;
  int counter_comp_mask=0;

  for(int j=0;j<matrix->getNumberOfColumns();j++){
    for(int i=0;i<matrix->getNumberOfRows();i++){
      if(comp_mask[counter_comp_mask]==1){
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


/* 
 * helper function for greedy_connectivity_1 sets certain rows and columns to -DBL_MAX 
 * @pram: pointer to the row which should be invalidated and turned to -inf 
 * @pram: pointer to the column which should be invalidated and turned to -inf
 * @pram: adjacency matrix for graph1
 * @pram: adjacency matrix for graph2
 * @pram: matrix indicating scores for nodal pairs 
 */
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


vector<int>* intersect(int* cols, int cols_size, struct coordinate_pair** row_cols,int row_cols_size){
  vector<int>* ret_value=new vector<int>();
  int counter=0;

  for(int j=0;j<row_cols_size;j++){ 
    int curr_row=row_cols[j]->row; 
    for(int i=0;i<cols_size;i++){
      if(curr_row==cols[i]){
	ret_value->insert(ret_value->begin()+counter,curr_row);
	counter++;
      }
    }
  }
  
  return ret_value;
}



int* get_valid_entries(SparseMatrix<float> graph1, int* ass,int size,int* ret_size){
  
  vector<int> assigned;
  SparseMatrix<float>* graph_copy=new SparseMatrix<float>(graph1);

  for(int i=0;i<size;i++){
    if(ass[i]!=-1){
      assigned.insert(assigned.begin()+i,i);
    }
  }
  int invalidate=1;
  *ret_size=0;

  for(int i=0;i<graph_copy->getNumberOfRows();i++){
    for(int k=0;k<assigned.size();k++){
      if(assigned[k]==i){
	invalidate=0;
	break;
      }
    }
    if(invalidate==1){
      for(int j=0;j<graph_copy->getNumberOfColumns();j++){
	(*graph_copy)[i][j]=-DBL_MAX;
      }
    }
    else{
      for(int j=0;j<graph_copy->getNumberOfColumns();j++){
	if((*graph_copy)[i][j]==1)
	  *ret_size=*ret_size+1;
      }
      invalidate=1;
    }
  }

  int* ret_arr=new int[*ret_size];
  int ret_arr_counter=0;
  int counter=0;


  for(int i=0;i<graph1.getNumberOfColumns();i++){
    for(int j=0;j<graph1.getNumberOfRows();j++){
      if((*graph_copy)[j][i]==1){
	ret_arr[ret_arr_counter]=i;
	ret_arr_counter++;
      }
      counter++;
    }
  }

  return ret_arr;

}


template <typename DT>
struct coordinate_pair** find_all_values(SparseMatrix<DT>& local_matches,DT* val,int val_size,int* row_cols_size){

  int rows=local_matches.getNumberOfRows();
  int cols=local_matches.getNumberOfColumns();

  int size=0;  
   SparseMatrix<DT>* local_matches_copy=new SparseMatrix<DT>(local_matches);

for(int i=0;i<local_matches_copy->getNumberOfRows();i++){
    for(int j=0;j<local_matches_copy->getNumberOfColumns();j++){
      for(int k=0;k<val_size;k++){
	if(val[k]==(*local_matches_copy)[i][j]){
	  (*local_matches_copy)[i][j]=-DBL_MAX;
	  size++;
	}
      }
    }
    }

    local_matches_copy=new SparseMatrix<DT>(local_matches);

    *row_cols_size=size;
 struct coordinate_pair **ret_value= new struct coordinate_pair*[size];
 struct coordinate_pair *pair;
 int counter=0,val_counter=0,ret_val_counter=0;

 for(int i=0;i<local_matches.getNumberOfRows();i++){
    for(int j=0;j<local_matches.getNumberOfColumns();j++) {
      if(val[val_counter]==counter) {
	pair=(struct coordinate_pair*)malloc(sizeof(struct coordinate_pair));
	pair->row=i;
	pair->col=j;
	ret_value[ret_val_counter]=pair;
	ret_val_counter++;
	val_counter++;
      }
      counter++;
    }
 }
 return ret_value;
} 






/*
 * checks whether all entries in a matrix are negative
 * @pram: an instance of a SparseMatrix that has all negative numbers
 */

template <typename DT>
int all_inf(SparseMatrix<DT> &mat){
  for(int i=0;i<mat.getNumberOfRows();i++){
    for(int j=0;j<mat.getNumberOfColumns();j++){
      if(mat[i][j]>0){
	return 0;
      }
    }
  }
  return 1;
}

/*
 * sets all entries in a matrix to be -inf
 * @pram: matrix that gets set to -inf
 */
template <typename DT>
void set_to_min(SparseMatrix<DT>& matrix){
 
  for(int i=0;i<matrix.getNumberOfRows();i++){
    for(int j=0;j<matrix.getNumberOfColumns();j++){
      matrix[i][j]=-DBL_MAX;
    }
  }

}

/*
 * sets certain values of matrix1 to be certain values of matrix2
 * @pram: SparseMatrix which gets changed
 * @pram: SparseMatrix whose values are used to change matrix1
 * @pram: vector representing the rows that need to be changed
 * @pram: vector representing the columns that need to be changed 
 */
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

/*
 * finds all the occurences >= a certain value in the sparse matrix
 * @pram: Sparse Matrix we read in
 * @pram: value that we are comparing
 * @pram: the number of times the value shows up
 */
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
  int retarr_counter=0;
  
  for(int i=0;i<matches2.getNumberOfRows();i++){
    for(int j=0;j<matches2.getNumberOfColumns();j++){
      if(matches2[i][j]>=value){
	 retarr[counter]=retarr_counter;
	  counter++;
	}
      retarr_counter++;
    }
    }
 
  return retarr;
} 

/*
 * returns the idth occurence of a value greater than val in the matrix
 * @pram: matrix that is being read in
 * @pram: the occurence of a value >= val we wish to return
 * @pram: pointer to row integer which gets set to the value we return
 * @pram: pointer to column integer which gets set to the value we return
 */
template <typename DT>
DT get_Max(SparseMatrix<DT> *matches, int id, DT val,int *row, int *col){
  int count=0;
  if(val<0){
    return -1;
  }


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


/*
 * returns the column of the largest value in a matrix, breaks ties randomly
 * @pram: matrix that represents nodal pair scores
 * @pram: pointer to variable that is used to indicate how good matching is
 * @pram: pointer to row variable which we set to the row of the largest value
 * @pram: pointer to column variable which we set to the column of the largest value
 */

template <typename DT>
int return_max(SparseMatrix<DT>& matches, DT* total_score,int* max_row,int* max_col){
  DT max_so_far=-DBL_MAX;
  int max_so_far_count=1;
  *max_row=0; 
  *max_col=0;

  for(int j=0;j<matches.getNumberOfColumns();j++){
    for(int i=0; i<matches.getNumberOfRows();i++){
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

  if(max_so_far==-DBL_MAX){
    return -1;
  }


  *total_score+=max_so_far;

  return *max_col;
}

/*
 * turns all values located in row row and column col into -inf
 * @pram: row we wish to turn into -inf
 * @pram: column we wish to turn into -inf
 * @pram: matrix on which we perform operation
 */
template <typename DT>
void invalidate(int row, int col, SparseMatrix<DT>& matches){
  for(int j=0;j<matches.getNumberOfColumns();j++){
    matches[row][j]=-DBL_MAX;
  }
  for(int j=0;j<matches.getNumberOfRows();j++){
    matches[j][col]=-DBL_MAX;
  }

}
/*
 * function called if matching is not complete (some nodes in graph1
 * don't get matched to any nodes in graph2)
 * @pram: an array which signifies the matching between graph1 and graph2
 * @pram: SparseMatrix representing graph1
 * @pram: SparseMatrix representing graph2 
 */
void match_rest(int* assignment, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2){
  printf("match rest called\n");

  if(graph1.getNumberOfRows()<=graph2.getNumberOfRows()){
    int unassigned_graph1[graph1.getNumberOfRows()];
    int unassigned_graph2[graph2.getNumberOfRows()];
    int counter1=0;
    for(int j=0;j<graph1.getNumberOfRows();j++){
      unassigned_graph1[j]=-1;
    }
    for(int j=0;j<graph2.getNumberOfRows();j++){
      unassigned_graph2[j]=1;
    }
    
    for(int j=0;j<graph1.getNumberOfRows();j++){
      if(assignment[j]==-1){
	unassigned_graph1[counter1]=j;
	counter1++;
      }
      else{
	unassigned_graph2[assignment[j]]=-1;
      }
    }
   
    for(int i=0;unassigned_graph1[i]!=-1;i++){
      for(int j=0;j<graph2.getNumberOfRows();j++){
	if(unassigned_graph2[j]==1){
	  assignment[unassigned_graph1[i]]=j;
	  unassigned_graph2[j]=-1;
	  break;
	}
      }
    }
    
  }
  else {
    int counter=graph2.getNumberOfRows();
    for(int i=0;i<graph1.getNumberOfRows();i++){
      if(assignment[i]==-1){
	assignment[i]=counter;
	counter++;
      }
    }
  }

}

template <typename DT>
SparseMatrix<DT>* getPermMatrix(const int size){
  int* assignment=(int*)malloc(sizeof(int)*size);
  SparseMatrix<DT>* ret_matrix=new SparseMatrix<DT>(size,size);
  
  for(int i=0;i<size;i++){
    ret_matrix[i][assignment[i]]=1;
  }
  return ret_matrix;
}



#endif

