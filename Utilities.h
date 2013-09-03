/*************************************************************************************************
 * This file contains functions that are used to perform transformations on matrices and vectors.*
 *                                                                                               *
 *************************************************************************************************/

#ifndef _UTILITIES_h
#define _UTILITIES_h

#include "Vertex.h"
#include "Matrices/DenseMatrix1D.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stack>


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
 * @pram: pointer to array
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

/*
 * returns the maximum value in the std::vector vec
 * @pram: pointer to std::vector
 */
template<typename DT>
std::vector<DT>* vector_max(std::vector<DT>* vec){
    DT max_so_far= (*vec)[0];
    std::vector<int>* ret_vec= new std::vector<int>();
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

/*
 * returns an array of integers which indicate
 * the indices where the value val is located in
 * the array
 * @param: array we are traversing through
 * @param: the value we are looking for
 * @param: the size of the array
 */
template <typename DT>
int* find_in_arr(DT* arr,DT val, int arr_size){
    int counter=0;
    
    //count the number of times val shows up in arr
    for(int i=0;i<arr_size;i++){
        if(arr[i]==val)
            counter++;
    }
    
    //create and return the array with indices of where
    //value is located
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
 * multiples a std::vector by a scalar
 * @pram: pointer to the array returned
 * @pram: pointer to array used to represent the std::vector
 * @pram: the number by which the std::vector gets scaled
 */
template <typename DT>
DT* scalar_multiplication(DT *old_row,int size, DT scaling_factor){
    
    DT* ret_row= new DT[size];
    
    //multiply each value in array old_row with value scaling_factor
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
 * @param: std::vector of pointers to vertex objects
 * @param: the component we wish to mask
 */
std::vector<int>* component_mask(std::vector<vertex*>& vertices, int component){
    
    std::vector<int>* comp_mask = NULL;
    bool hasConnectedComponeents = false;
    
    for(int i=0;i< vertices.size();i++)
    {
        vertex* curr_vertex= vertices[i];
        if(curr_vertex->get_low_link()==component){
            if (!hasConnectedComponeents)
            {
                comp_mask = new std::vector<int>(vertices.size(),0);
                hasConnectedComponeents = true;
            }
            (*comp_mask)[i] = 1;
        }
    }
    return comp_mask;
}

#endif
