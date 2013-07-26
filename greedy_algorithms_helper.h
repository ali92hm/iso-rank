#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <cmath>
#include "Matricies/DenseMatrix.h"
#include <limits>

#ifndef greedy_algorithm_helper_h
#define greedy_algorithm_helper_h


vector<int>* intersect(int*, int, struct coordinate_pair**,int);
void match_rest(int*, DenseMatrix<float>&, DenseMatrix<float>&);
int* get_valid_entries(DenseMatrix<float>, int*,int,int*);
vector<int>* choose_cols(struct coordinate_pair**,int,int);

struct coordinate_pair{
  int row;
  int col;
};



/*
 * compares two floats and returns whether 
 * an integer to indicate which is bigger
 * 0 if a=b, 1 if a>b, -1 a<b
 * @pram: float to be compared
 * @pram: float to be compared
 */
int compareFloats(float a, float b){
  float smallest_float=numeric_limits<float>::epsilon();
  
  if(abs(a-b)<smallest_float){
    return 0;
  }
  else if(a-b>0)
    return 1;
  else
    return -1;

}


/*
 * returns a vector of the columns col 
 * such that (row,col) is a pair in rc
 * @pram:array of structs of coordinate pairs
 * @pram: size of the array of structs
 * @pram: the row we wish to match
 */
vector<int>* choose_cols(struct coordinate_pair** rc,int rc_size,int row){
  vector<int> *ret_vec= new vector<int>(0);

  for(int i=0;i<rc_size;i++){
    struct coordinate_pair *pair=rc[i];
    
    if(pair->row==row){
      ret_vec->push_back(pair->col);
    }
  }
  
  return ret_vec;
}


/*
 * returns a DenseMatrix Object which is a reshaped eigenvector
 * @pram: a pointer to an array of doubles which represents the eigenvector
 * @pram: number of rows in the matrix returned
 * @pram: number of columns in the matrix returned
 * @pram: component mask vector indicating which nodes are present in current component
 */
template <typename DT>
DenseMatrix<DT> reshape(DT* eigenvector,const int rows,const int cols, vector<int> &comp_mask){

  DenseMatrix<DT> matrix(rows,cols);
  int counter_eig_vector=0;
  int counter_comp_mask=0;


    for(int i=0;i<matrix.getNumberOfRows();i++){
        for(int j=0;j<matrix.getNumberOfColumns();j++){
      if(comp_mask[counter_comp_mask]==1){
	matrix[i][j]=eigenvector[counter_eig_vector];
	counter_eig_vector++;
      }
      else{
	matrix[i][j]=0;
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
void neighbor_enforcement(int* row_index,int* col_index, DenseMatrix<float>& graph1,DenseMatrix<float>& graph2, DenseMatrix<DT>& matches){

  
   for(int i=0;i<graph1.getNumberOfColumns();i++){
      if(graph1[*row_index][i]==1){
	for(int j=0;j<graph2.getNumberOfRows();j++){
	  //if node i neighbors node row_index in graph1
	  // and node j does not neighbor col_index invalidate (i,j) matching
	  if(graph2[j][*col_index]==0){
	    matches[i][j]=-DBL_MAX;
	  }
	}
      }
    }
  
    for(int i=0;i<graph2.getNumberOfColumns();i++){
      if(graph2[*col_index][i]==1){
	for(int j=0;j<graph1.getNumberOfRows();j++){
	    //if node j neighbors node col_index in graph2
	  // and node i does not neighbor row_index invalidate (i,j) matching
	  if(graph1[j][*row_index]==0){
	    matches[j][i]==-DBL_MAX;
	  }
	}
      }
    }
}

/*
 * returns a vector of nodes from graph1 such that the 
 * node r exists as a pair (r,c) in row_cols and 
 * c is exists in array cols
 * @param: array of nodes in graph2 being considered for the matching
 * @param: size of array that contains nodes from graph2
 * @param: array of coordinate_pair structs being considered for matching
 * @param: size of the array row_cols
 */
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

/*
 * returns array of nodes from graph2 that can be made availabe for matching
 * @param: DenseMatrix representing graph1
 * @param: array of assignments
 * @param: size of assignments array
 * @param: size of returned array
 */
int* get_valid_entries(DenseMatrix<float> graph1, int* ass,int size,int* ret_size){
  
  vector<int> assigned;
  DenseMatrix<float>* graph_copy=new DenseMatrix<float>(graph1);

  for(int i=0;i<size;i++){
    if(ass[i]!=-1){
      assigned.push_back(i);
    }
  }
  int invalidate=1;
  *ret_size=0;

  //invalidate certain nodes from consideration for matching
  for(int i=0;i<graph_copy->getNumberOfRows();i++) {
   
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


    for(int j=0;j<graph1.getNumberOfRows();j++){
      for(int i=0;i<graph1.getNumberOfColumns();i++){
      if((*graph_copy)[j][i]==1){
	ret_arr[ret_arr_counter]=i;
	ret_arr_counter++;
      }
      counter++;
    }
  }
  delete graph_copy;

  return ret_arr;

}


/*
 * returns an array of coordinate_pair structs (r,c)
 * such that r,c in scores_matrix has a high score
 * @param: DenseMatrix that represents the nodal pairings scores
 * @param: array of values that are high enough and in local_matches
 * @param: size of array val
 * @param: size of the array returned
 */
template <typename DT>
struct coordinate_pair** find_all_values(DenseMatrix<DT>& local_matches,DT* val,int val_size,int* row_cols_size){

  int rows=local_matches.getNumberOfRows();
  int cols=local_matches.getNumberOfColumns();

  int size=0;  
   DenseMatrix<DT>* local_matches_copy=new DenseMatrix<DT>(local_matches);

   //count the number of values that we wish to consider
   for(int i=0;i<local_matches_copy->getNumberOfRows();i++){
    for(int j=0;j<local_matches_copy->getNumberOfColumns();j++){
      for(int k=0;k<val_size;k++){
	if(compareFloats(val[k],(*local_matches_copy)[i][j])==0){
	  (*local_matches_copy)[i][j]=-DBL_MAX;
	  size++;
	}
      }
    }
   }

 delete local_matches_copy;

    local_matches_copy=new DenseMatrix<DT>(local_matches);

    *row_cols_size=val_size;
 struct coordinate_pair **ret_value= new struct coordinate_pair*[val_size];
 struct coordinate_pair *pair;
 int counter=0,val_counter=0,ret_val_counter=0;

 //add all the nodal pairings that are greater than a certain value to the final array of pairings
 for(int i=0;i<local_matches.getNumberOfRows();i++){
    for(int j=0;j<local_matches.getNumberOfColumns();j++) {
      if(val_counter<val_size&&val[val_counter]==counter) {
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
 delete local_matches_copy;
 return ret_value;
} 






/*
 * checks whether all entries in a matrix are negative
 * @pram: an instance of a DenseMatrix that has all negative numbers
 */

template <typename DT>
int all_inf(DenseMatrix<DT> &mat){
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
void set_to_min(DenseMatrix<DT>& matrix){
 
  for(int i=0;i<matrix.getNumberOfRows();i++){
    for(int j=0;j<matrix.getNumberOfColumns();j++){
      matrix[i][j]=-DBL_MAX;
    }
  }

}

/*
 * sets certain values of matrix1 to be certain values of matrix2
 * @pram: DenseMatrix which gets changed
 * @pram: DenseMatrix whose values are used to change matrix1
 * @pram: vector representing the rows that need to be changed
 * @pram: vector representing the columns that need to be changed 
 */
template<typename DT>
void set_matrix_values(DenseMatrix<DT>& matrix1,DenseMatrix<DT>& matrix2,  vector<int>& rows, vector<int>& cols){
  int r;
  int c;
  

  for(int i=0;i<matrix1.getNumberOfRows();i++){
    for(int j=0;j<matrix1.getNumberOfColumns();j++){
      matrix1[i][j]=-DBL_MAX;
    }
  }


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
DT* find_values(DenseMatrix<DT>& matches2,DT value,int *size){
  *size=0;

  //increments a counter to count the number of values in the matrix
  //greater than a certain amount
  for(int i=0; i<matches2.getNumberOfRows();i++){
    for(int j=0;j<matches2.getNumberOfColumns();j++){
      if(compareFloats(matches2[i][j],value)==1||
	 compareFloats(matches2[i][j],value)==0) {
	*size=*size+1;
      }
    }
  }

  DT* retarr= new DT[*size];
  
  int counter=0;
  int retarr_counter=0;
  
  //returns an array indicating where in the matrix those values 
  //exist
  for(int i=0;i<matches2.getNumberOfRows();i++){
    for(int j=0;j<matches2.getNumberOfColumns();j++){
      if(compareFloats(matches2[i][j],value)==1||
	 compareFloats(matches2[i][j],value)==0){
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
DT get_Max(DenseMatrix<DT> *matches, int id, DT val,int *row, int *col){
  int count=0;
  if(val<0){
    return -1;
  }

  for(int i=0; i<matches->getNumberOfRows();i++){
    for(int j=0;j<matches->getNumberOfColumns();j++){
      if(compareFloats((*matches)[i][j],val)==1||compareFloats((*matches)[i][j],val)==0){
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
int return_max(DenseMatrix<DT>& matches, DT* total_score,int* max_row,int* max_col){

  DT max_so_far=-DBL_MAX;
  int max_so_far_count=1;
  *max_row=0; 
  *max_col=0;


  //keeps track of the largest value in the scores matrix
  //and keeps a count of the number of times the values occurs
    for(int i=0; i<matches.getNumberOfRows();i++){
      for(int j=0;j<matches.getNumberOfColumns();j++){
      if(compareFloats(matches[i][j],max_so_far)==0){
	max_so_far_count++;
      }
      else if(compareFloats(matches[i][j],max_so_far)==1){
	max_so_far=matches[i][j];
	max_so_far_count=1;
      }  
    }

  }
 
    //choose a random occurence of the maximum value
    // to set max_row and max_col to

  int counter=0;
  int random_number= rand()%max_so_far_count+1;
  int set=0;

  for(int i=0; i<matches.getNumberOfRows();i++){
    for(int j=0;j<matches.getNumberOfColumns();j++){
      if(compareFloats(matches[i][j],max_so_far)==0){
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

  if(compareFloats(max_so_far,-DBL_MAX)==0) {
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
void invalidate(int row, int col, DenseMatrix<DT>& matches){
  //set every entry in row and col to be -inf
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
 * @pram: DenseMatrix representing graph1
 * @pram: DenseMatrix representing graph2 
 */
void match_rest(int* assignment, DenseMatrix<float>& graph1, DenseMatrix<float>& graph2){
  
  if(graph1.getNumberOfRows()<=graph2.getNumberOfRows()){
    int unassigned_graph1[graph1.getNumberOfRows()];
    int unassigned_graph2[graph2.getNumberOfRows()];
    int counter1=0;

    //intialize arrays
   for(int j=0;j<graph1.getNumberOfRows();j++){
      unassigned_graph1[j]=-1;
    }
    for(int j=0;j<graph2.getNumberOfRows();j++){
      unassigned_graph2[j]=1;
    }
    
    //assign unassigned nodes
    for(int j=0;j<graph1.getNumberOfRows();j++){
      if(assignment[j]==-1){
	unassigned_graph1[counter1]=j;
	counter1++;
      }
      else{
	unassigned_graph2[assignment[j]]=-1;
      }
    }
    int node_counter=graph2.getNumberOfRows();
    for(int i=0;unassigned_graph1[i]!=-1;i++){
      int j=0;
      for(;j<graph2.getNumberOfRows();j++){
	if(unassigned_graph2[j]==1){
	  assignment[unassigned_graph1[i]]=j;
	  unassigned_graph2[j]=-1;
	  break;
	}
      }
    }


    for(int i=0;i<graph1.getNumberOfRows();i++){
      if(assignment[i]==-1){
	assignment[i]=node_counter;
	node_counter++;
      }
    }   
  }
  else {
    int unassigned_graph1[graph1.getNumberOfRows()];
    int unassigned_graph2[graph2.getNumberOfRows()];
    int counter1=0;

    //intialize arrays
   for(int j=0;j<graph1.getNumberOfRows();j++){
      unassigned_graph1[j]=-1;
    }
    for(int j=0;j<graph2.getNumberOfRows();j++){
      unassigned_graph2[j]=1;
    }
    
    //assign unassigned nodes
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

    int counter=graph2.getNumberOfRows();
    for(int i=0;i<graph1.getNumberOfRows();i++){
      if(assignment[i]==-1){
		assignment[i]=counter;
		counter++;
      }
    }
  }

}


/*
 * returns a permutation matrix with dimensions size x size
 * @param: size of the permutation matrix
 */
DenseMatrix<float> getPermMatrix(int *ass, int ass_size,const int bigger_matrix_size){
  DenseMatrix<float> ret_matrix(bigger_matrix_size,bigger_matrix_size);
  int hold;
  for(int i=0;i<bigger_matrix_size;i++){
    if(i<ass_size){
       hold=ass[i];
       ret_matrix[i][hold]=1;
    }
    else{
      ret_matrix[i][i]=1;
    }
    
  }
  return ret_matrix;
}

/*
 * looks through the vector of nodes and removes the 
 * ones that have already been assigned
 * @param: array showing the nodal assignments
 * @param: vector filled with nodes that are neighbors of a node 
 */
void invalidate_neighbors(int* assignment, vector<int> neigh)
{
   int hold;

   
   for(int i=0;i<neigh.size();i++)
   {
      hold=neigh[i];
      if(assignment[hold]==1)
      {
		neigh.erase(neigh.begin()+i);
      }
    }

}

/*
 *initializes an array to have all indices set to init_val
 *@param: array we wish to initialize
 *@param: size of the array we wish to initiliaze
 *@param: initial value that array should be populated with
 */
template <typename DT>
void init_array(DT* arr, int arr_size,DT init_val){
  
  for(int i=0;i<arr_size;i++){
    arr[i]=init_val;
  }
  

}



#endif
