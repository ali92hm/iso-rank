#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <cmath>
#include "SparseMatrix.h"
#include <limits>
#ifndef greedy_algorithm_h
#define greedy_algorithm_h



struct coordinate_pair{
  int row;
  int col;
};

vector<int>* intersect(int*, int, struct coordinate_pair**,int);
void match_rest(int*, SparseMatrix<float>&, SparseMatrix<float>&);
int* get_valid_entries(SparseMatrix<float>, int*,int,int*);
vector<int>* choose_cols(struct coordinate_pair**,int,int);

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

  //intialize assignment array
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
  }

  for(int i=0;i<min(graph1_nodes,graph2_nodes);i++){

    //get maximum score in matrix and set assignment
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

  //intialize assignment array 
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
  }
 
  
  for(int i=0;i<min(graph1_nodes,graph2_nodes);i++){
    
    //find maximum in scores matrix and perform assignment 
    return_max(matches,&total_score,&row,&col);
    assignment[row]=col;
    invalidate(row,col,matches);
   
    //change matrix s.t. only neighbors of row are allowed to 
    //match to neighbors of col
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

  //intialize assignment array
  for(int i=0;i<graph1_nodes;i++){
    assignment[i]=-1;
    assignment2[i]=0;
  }

  SparseMatrix<DT>* active_matches=new SparseMatrix<DT>(matches);
  DT* idxarr;
  vector<int> assigned_G1;
  int size=0,random_id,vector_size=0,curr_row;
 
srand(time(NULL));

//run while loop until all nodes are assigned and scores matrix isn't all negative
 while(sum_array(assignment2,graph1_nodes)<min(graph1_nodes,graph2_nodes)&& return_max(*active_matches,&score,&row,&col)>-1)
{

      if(all_inf(*active_matches)){
        match_rest(assignment,graph1,graph2); 
        printf("score of matching: %f\n", final_score);
	return;
      }
 
      //find all values in scores matrix greater than a certain amount 
      idxarr=find_values(*active_matches, score-prev_score-max_tol,&size);
      random_id=rand()%size+1;
      max= get_Max(active_matches,random_id,score-prev_score-max_tol,&row,&col);
    
      //perform assignment by choosing a random pair thats high enough
      final_score+=max;
      assignment[row]=col;
      assignment2[row]=1;
      
      assigned_G1.push_back(row);
      vector_size++;
      invalidate(row, col, matches);
      set_to_min(*active_matches);
      int curr_col;

      //change scores matrix such that only neighbors of already matched
      //nodes are allowed to match to one another
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
  
  //intialize all arrays 
  for(int i=0;i<graph1_nodes;i++) {
    assignment[i]=-1;
    assignment2[i]=0;
    assignment_G1[i]=0;
    assignment_G2[i]=0;
  }

  vector<int>* neigh_1;
  vector<int>* neigh_2;

  //run while loop until all nodes are assigned 
  while(sum_array(assignment2,graph1_nodes)<min(graph1_nodes,graph2_nodes)){
  
    //find the highest matching score and make that assignment 
    return_max(matches, &final_score,&row,&col);
    assignment[row]=col;
    assignment2[row]=1;
    assignment_G1[row]=1;
    assignment_G2[col]=1;
    invalidate(row,col,matches);



    //find neighbors of already assigned nodes and invalidate already assigned nodes from further consideration
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

    //if scores matrix is all -inf match unassigned nodes and return 
     if(all_inf(matches)){
       match_rest(assignment,graph1,graph2);
       return;
     }

     //run for loop until all neighbors are assigned and score matrix isn't all -inf
    for(int i=0;i<min(neigh_1->size(),neigh_2->size())&&!all_inf(*local_matches);i++){
      
      //find best nodal pairing and peform assignment 
      return_max(*local_matches,&final_score,&row,&col);
      assignment[row]=col;
      assignment2[row]=1;
      assignment_G1[row]=1;
      assignment_G2[col]=1;

      //invalidate already assigned nodes from further consideration
      invalidate(row,col,*local_matches);
      invalidate(row,col,matches);

       //if scores matrix is all -inf match unassigned nodes and return 
        if(all_inf(matches)){
	  match_rest(assignment,graph1,graph2);
	  return;
	}

    }
  }
  return;
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
 * performs a greedy matching and enforces connectivity by proceeding outwards radially
 * chooses the most connected neighbor at every iteration
 * @pram:matrix indicating scores for nodal pairs 
 * @pram: adjacency matrix for graph1
 * @pram: adjacency matrix for graph2
 * @pram: pointer to the array that indicates the best matching 
 */
template <typename DT>
void greedy_connectivity_4(SparseMatrix<DT>& matches, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2,int* assignment){
  
DT final_score=0;
int* add_order=new int[graph1.getNumberOfRows()];
int* ass=new int[graph1.getNumberOfRows()];
int add_idx=0;
DT score=0;
DT max_tol=pow(10,-6);
int row,col,size=0;
int assigned_G1[graph1.getNumberOfRows()];
int assigned_G2[graph2.getNumberOfRows()];
vector<int>* neigh_1;
vector<int>* neigh_2;
int add_order_counter=2;


 //intializing all arrays 
 for(int j=0;j<graph1.getNumberOfRows();j++){
   add_order[j]=-1;
   ass[j]=0;
   assigned_G1[j]=-1;
 }

 for(int j=0;j<graph2.getNumberOfRows();j++){
      assigned_G2[j]=-1;
 }

  
 return_max(matches, &score,&row,&col);

DT* idx_array =find_values(matches,score - max_tol,&size);
int random_id=rand()%size+1; 
DT max= get_Max(&matches,random_id,score-max_tol,&row,&col);
 
 //assign first row column pair
 final_score+=max;
 assignment[row]=col; 
 assigned_G1[row]=1;
 assigned_G2[col]=1;
 invalidate(row,col,matches);
 add_order[0]=row;


 neigh_1= graph1.getNeighbors(row);
 neigh_2= graph2.getNeighbors(col);
 delete []idx_array;


 int hold;

 //remove nodes that are already assigned from consideration in the next assignment
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


 //create local score matrix and only set the values for the pairs we're considering 
 SparseMatrix<DT>* matches_local=new SparseMatrix<DT>(matches.getNumberOfRows(),matches.getNumberOfColumns());
 set_matrix_values(*matches_local, matches, *neigh_1, *neigh_2);

 score=0;
 size=0;
 return_max(*matches_local, &score,&row,&col);
 idx_array=find_values(*matches_local,score-max_tol,&size);
 random_id=rand()%size+1;
 max=get_Max(matches_local,random_id,score-max_tol,&row,&col);

 //assign the second row and column pair
 final_score+=max;
 assignment[row]=col;
 assigned_G1[row]=1;
 assigned_G2[col]=1;
 invalidate(row,col,*matches_local);
 invalidate(row,col,matches);
 add_order[1]=row;

 delete []idx_array;

int counter=0;

 //remove nodes that are already assigned from consideration in the next assignment
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

    int best_row=0,best_col=0;
    
    //while loop that runs until either the last node is assigned or we run out of possible matchings
    while(add_order[graph1.getNumberOfRows()-1]==-1) {
    
      //for loop that aims to match all the neighbors of the currently selected nodal pairing
      for(int s=0; s<min(neigh_1->size(),neigh_2->size());s++) {
       	score=0;
       
	//finds all node pairings that are above a certain score and stores them in array rows_cols
	return_max(*matches_local,&score,&row,&col);
	idx_array=find_values(*matches_local,score-max_tol,&size);
	int rows_cols_size;
	int valid_entries_size=0,valid_entries2_size=0;
	struct coordinate_pair **rows_cols=find_all_values(*matches_local,idx_array,size,&rows_cols_size);
	delete []idx_array;
 
	//if number of nodal pairings with a high score is greater than 1
     	if(size>1) {

	  //looks through graph1 to see which edges exist and their intersection with node pairings with high enough scores
	  int* valid_entries= get_valid_entries(graph1,assignment,graph1.getNumberOfRows(),&valid_entries_size);
	  vector<int>* prev_assigned = intersect(valid_entries,valid_entries_size,rows_cols,size);
	  vector<int> g1c_count(graph1.getNumberOfRows());
	  int g1c_count_counter=0;

	  //find the connectivity of nodes being considered for matching
	  for(int f=0;f<graph1.getNumberOfRows();f++){
	    g1c_count[f]=-1;
	  }

	  for(int k=0;k<prev_assigned->size();k++){
	    int id=(*prev_assigned)[k];
	    int sum=0;

	    for(int j=0;j<valid_entries_size;j++){
	      if(valid_entries[j]==id) {
		sum++;
	      }
	    }

	    g1c_count[id]=sum;
	    g1c_count_counter++;
	  }
	  
	  //find the node best_row with highest connectivity to match
	  vector<int> *max_g1c=vector_max(&g1c_count);
	  int rand_number = rand()%(max_g1c->size());
	  best_row=(*max_g1c)[rand_number];
	  vector<int>* best_cols= new vector<int>();
	  int best_cols_counter=0;

	  //find all nodes in graph2 that are available for matching to 
	  //node just chosen from graph1
	  for(int r=0;r<rows_cols_size;r++) {
	    struct coordinate_pair *c_pair=rows_cols[r];
	    if(c_pair->row==best_row) {
	      best_cols->push_back(c_pair->col);
	      best_cols_counter++;
	    }
	  }


	  //find node from graph2 to match to node just chosen from graph1 
	  int* valid_entries2= get_valid_entries(graph2,assignment,graph1.getNumberOfRows(),&valid_entries2_size);
	  vector<int> *cols_chosen=choose_cols(rows_cols,rows_cols_size,best_row);

	  vector<int> g2c_count(graph2.getNumberOfRows());
	  int g2c_count_counter=0;

	  //finds the connectivity of each of the nodes being considered
	  for(int k=0;k<cols_chosen->size();k++) {
	    int id=(*cols_chosen)[k];
	    int sum=0;

	    for(int m=0;m<valid_entries2_size;m++){
	      if(valid_entries2[m]==id){
		sum++;
	      }
	    }
	    g2c_count[id]=sum;
	    g2c_count_counter++;
	  }

	  vector<int> *max_g2c=vector_max(&g2c_count);
	   rand_number = rand()%(max_g2c->size());
	   best_col=(*max_g2c)[rand_number];
	   
	     
	delete best_cols;
	delete max_g1c;
	delete max_g2c;
	delete prev_assigned;
	delete []valid_entries2;
	delete []valid_entries;
	delete cols_chosen;

	}
	else if(size==1){
	  //if number of pairings is just 1
	  struct coordinate_pair *cp=rows_cols[0];
	  best_row=cp->row;
	  best_col=cp->col;
     	}
	else{	 
	  for(int x=0;x<rows_cols_size;x++){
	    free(rows_cols[x]);
	  }
      	
	  delete []rows_cols; 
	  break;
	}

	//perform the assignment and invalidate corresp. rows and columns in scores matrix 
	//from further consideration
       	assignment[best_row]=best_col;
	assigned_G1[best_row]=1;
	assigned_G2[best_col]=1;

	DT max_matches=matches[best_row][best_col];
	add_order[add_order_counter]=best_row;
	add_order_counter++;
	invalidate(best_row,best_col,*matches_local);
	invalidate(best_row,best_col,matches);
	 
	
	final_score+=max_matches;
	

	for(int x=0;x<rows_cols_size;x++){
	  free(rows_cols[x]);
	}

      	delete []rows_cols; 

	} //for min(neigh1, neigh2) 

      
	add_idx++;	
        int r=add_order[add_idx];

	//if a match not made at add_ixth iteration break
	if(r==-1){
	  add_order[graph1.getNumberOfRows()-1]=-2;
	  break;
	}

	//choose the next set of nodes to match   
	int c=assignment[r];
	delete neigh_1;
	delete neigh_2;

      	neigh_1=graph1.getNeighbors(r);
	neigh_2=graph2.getNeighbors(c);
	

	//remove nodes that have already been assigned from consideration
	 for(int a=0;a<neigh_1->size();a++){
	   hold=(*neigh_1)[a];
	   if(assigned_G1[hold]==1){
	     neigh_1->erase(neigh_1->begin()+a);
	     a--;
	   }
	 }

	 for(int b=0;b<neigh_2->size();b++){
	   hold=(*neigh_2)[b];
	   if(assigned_G2[hold]==1){
	     neigh_2->erase(neigh_2->begin()+b);
	     b--;
	   }
	 }

    
 
     set_matrix_values(*matches_local, matches, *neigh_1, *neigh_2);
     counter++;


      } //while add_order
    
    //if matching is not complete match the rest of the nodes
    for(int i=0;i<graph1.getNumberOfRows();i++){
      if(assignment[i]==-1){
	  match_rest(assignment,graph1,graph2);
      break;
    }


   }
    delete neigh_1;
    delete neigh_2;
    delete matches_local;
    delete []ass;
    delete []add_order;
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


    for(int i=0;i<matrix->getNumberOfRows();i++){
        for(int j=0;j<matrix->getNumberOfColumns();j++){
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
 * @param: SparseMatrix representing graph1
 * @param: array of assignments
 * @param: size of assignments array
 * @param: size of returned array
 */
int* get_valid_entries(SparseMatrix<float> graph1, int* ass,int size,int* ret_size){
  
  vector<int> assigned;
  SparseMatrix<float>* graph_copy=new SparseMatrix<float>(graph1);

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
 * @param: SparseMatrix that represents the nodal pairings scores
 * @param: array of values that are high enough and in local_matches
 * @param: size of array val
 * @param: size of the array returned
 */
template <typename DT>
struct coordinate_pair** find_all_values(SparseMatrix<DT>& local_matches,DT* val,int val_size,int* row_cols_size){

  int rows=local_matches.getNumberOfRows();
  int cols=local_matches.getNumberOfColumns();

  int size=0;  
   SparseMatrix<DT>* local_matches_copy=new SparseMatrix<DT>(local_matches);

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

    local_matches_copy=new SparseMatrix<DT>(local_matches);

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
DT* find_values(SparseMatrix<DT>& matches2,DT value,int *size){
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
DT get_Max(SparseMatrix<DT> *matches, int id, DT val,int *row, int *col){
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
int return_max(SparseMatrix<DT>& matches, DT* total_score,int* max_row,int* max_col){

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

srand (time(NULL));

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
void invalidate(int row, int col, SparseMatrix<DT>& matches){
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
 * @pram: SparseMatrix representing graph1
 * @pram: SparseMatrix representing graph2 
 */
void match_rest(int* assignment, SparseMatrix<float>& graph1, SparseMatrix<float>& graph2){

  
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


/*
 * returns a permutation matrix with dimensions size x size
 * @param: size of the permutation matrix
 */
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

