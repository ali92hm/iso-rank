//
//  util.h
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/14/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef _util_h
#define _util_h

#include "vertex.h"
#include "Matrix2D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stack>


using namespace std;


/*
 * returns the sum of the values in array arr
 * @pram: pointer to an array of floats
 * @pram: size of the array
 */
template <typename DT>
DT sum_array(DT* arr,int arr_size){
    DT sum=0;
    
    for(int j=0;j<arr_size;j++){
        sum+=arr[j];
    }
    
    return sum;
}


/*
 * returns the minimum value in the array arr
 * @pram: pointer to array of floats
 * @pram: size of the array
 */
template <typename DT>
DT min(DT* arr, int arr_size){
    DT min_so_far=arr[0];
    
    for(int i=1;i<arr_size;i++){
        if(arr[i]<min_so_far){
            min_so_far=arr[i];
        }
    }
    
    return min_so_far;
}


/*
 * returns the maximum value in the array arr
 * @pram: pointer to array of floats
 * @pram: size of the array
 */
template <typename DT>
DT max(DT* arr, int arr_size){
    DT max_so_far=arr[0];
    
    for(int i=1;i<arr_size;i++){
        if(arr[i]>max_so_far){
            max_so_far=arr[i];
        }
    }
    
    return max_so_far;
}


template<typename DT>
vector<DT>* vector_max(vector<DT>* vec){
  DT max_so_far= (*vec)[0];
  vector<int>* ret_vec= new vector<int>();
  int vec_counter=0;

  for(int i=0;i<vec->size();i++){
    if((*vec)[i]>max_so_far)
      max_so_far=(*vec)[i];
  }

  for(int i=0;i<vec->size();i++){
    if((*vec)[i]==max_so_far){
      ret_vec->insert(ret_vec->begin()+vec_counter,i);
      vec_counter++;
    }
  }
  return ret_vec;
}


/*
 * returns the mean of all the floats in an array
 * @pram: pointer to array of floats
 * @pram: size of array
 */
template <typename DT>
DT mean(DT *arr,int arr_size){
    return sum_array(arr,arr_size)/arr_size;
}



/*
 * returns the standard deviation of the values in arr
 * @pram: pointer to array of floats
 * @pram: size of the array
 */
template <typename DT>
DT std_dev(DT *arr, int arr_size){
    DT arr_mean=mean(arr,arr_size);
    
    DT variance=0;
    for(int index=0;index<arr_size;index++){
        DT hold =arr_mean-arr[index];
        variance+= hold*hold;
    }
    variance=variance/arr_size;
    return sqrtf(variance);
    
}


template <typename DT>
int* find_in_arr(DT* arr,DT val, int arr_size){
  int counter=0;
  
  for(int i=0;i<arr_size;i++){
    if(arr[i]==val)
      counter++;
  }
  int* ret_arr= new int[counter];
  counter=0;
  for(int i=0;i<arr_size;i++){
    if(arr[i]==val){
      ret_arr[counter]=i;
      counter++;
    }
  }

  return ret_arr;
}

/*
 * multiples a vector by a scalar
 * @pram: pointer to the array returned
 * @pram: pointer to array used to represent the vector
 * @pram: the number by which the vector gets scaled
 */
template <typename DT>
DT* scalar_multiplication(DT *old_row,int size, DT scaling_factor){
    
    DT* ret_row= new DT[size];
    
#pragma omp for
    for(int i=0;i<size;i++){
        ret_row[i]=scaling_factor*old_row[i];
    }
    
    return ret_row;
}

/*
 * takes an array of vertex objects vertices and an integer
 * signifying a component of the graph. The function returns
 * an integer array where vertices that don't belong to the
 * component are masked out.
 */
vector<int>* component_mask(vector<vertex*>& vertices, int component){
    
    vector<int>* comp_mask = NULL;
    bool hasConnectedComponeents = false;
    
    for(int i=0;i< vertices.size();i++)
    {
        vertex* curr_vertex= vertices[i];
        if((*curr_vertex).get_low_link()==component){
            if (!hasConnectedComponeents)
            {
                comp_mask = new vector<int>(vertices.size(),0);
                hasConnectedComponeents = true;
            }
            (*comp_mask)[i]= 1;
        }
    }
    return comp_mask;
}



#endif
