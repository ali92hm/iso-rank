//
//  Tarjan.h
//  Sparse_Matrix
//
//  Created by Ali Hajimirza on 6/14/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef _Tarjan_h
#define _Tarjan_h

#include "Matricies/DenseMatrix.h"
#include "Matricies/SparseElement.h"
#include "vertex.h"
#include <stack>

/*
 * performs a binary search on the edges of sparse graph
 * to find edge with curr_row as an endpoint
 * @pram: pointer to the 2-d array rep of a sparsegraph
 * @pram: int signifying the vertex we wish to work with
 */

template <typename T>
int binary_search_index(std::vector<SparseElement<T> >& sparse_graph, int curr_row){
    int start=0;
    int sparse_size = sparse_graph.size();
    int finish=sparse_size-1;
    int mid=0;
    while(start<=finish && mid>=0 && mid<sparse_size){
        mid=(start+finish)/2;
        
        if(sparse_graph[mid].getI() == curr_row){
            break;
        }
        else if(sparse_graph[mid].getI() < curr_row){
            start=mid+1;
        }
        else{
            finish=mid-1;
        }
    }
    
    while(mid>-1&&sparse_graph[mid].getI()==curr_row){
        mid--;
    }
    mid++;
    
    if(start>finish)
        return -1;
    
    return mid;
    
}

/*
 * returns the smaller of two numbers
 * @pram: first integer
 * @pram: second integer
 */
int min(int a, int b){
    if(a>b)
        return b;
    else
        return a;
    
}


/*
 * returns an integer signifying whether the vertex w
 * exists in the stack s
 * @pram: stack s used to perform tarjan's
 * @pram: vertex object w
 */
int contains(std::stack<vertex*>& s,vertex& w)
{
    std::stack<vertex*> other_stack;
    int found=0;
    vertex* current;
    
    while(!s.empty()){
        current = s.top();
        if((*current)==w){
            found=1;
            break;
        }
        s.pop();
        other_stack.push(current);
    }
    
    while(!other_stack.empty()){
        current=other_stack.top();
        other_stack.pop();
        s.push(current);
    }
    
    return found;
}
/*
 * Perform's tarjan's scc algorithm on the sparse graph given
 * @pram: pointer to array of pointers of sparse_matrix_element structs
 * @pram: the number of sprase_matrix_element's in sparse_graph
 * @pram: the number of vertices in the graph
 * @pram: pointer to index which is used in tarjan's algorithm to keep track of scc
 * @pram: the current vertex that tarjan's is being performed on
 * @pram: pointer to array of compIds
 * @pram: pointer to vertex objects
 * @pram: stack used to perform the algorithm
 */

template <typename T>
void strong_com(std::vector<SparseElement<T> >& sparse_graph, int num_vertices,int *index,int vertex_number,std::vector<vertex*>& vertices,std::stack<vertex*>& st){
    //std::cout << "calling strong component on" << vertex_number << std::endl;
    
    
    int sparse_size=sparse_graph.size();
    int sparse_edges_index = binary_search_index(sparse_graph,vertex_number);
    vertex* curr_vertex=vertices[vertex_number];
    vertex* other_vertex;
    SparseElement<T> sparse_vertex;
    curr_vertex->set_index(*index);
    curr_vertex->set_low_link(*index);
    (*index)=(*index)+1;
    //cout << "index: " << *index << endl;
    
    if( (sparse_edges_index<sparse_size) && (sparse_edges_index > -1) ){
        sparse_vertex=sparse_graph[sparse_edges_index];
    }
    else{
        return;
    }
    
    
    st.push(curr_vertex);
    while( (sparse_edges_index < sparse_size) && (sparse_vertex.getI()==vertex_number) ){
	sparse_vertex=sparse_graph[sparse_edges_index];       
	 other_vertex=vertices[sparse_vertex.getJ()];
        
        if(other_vertex->get_index()==-1){
            strong_com(sparse_graph, num_vertices,index,sparse_vertex.getJ(),vertices,st);
            curr_vertex->set_low_link(min(curr_vertex->get_low_link(),other_vertex->get_low_link()));
        }
        else if(contains(st,(*other_vertex))==1){
            curr_vertex->set_low_link(min(curr_vertex->get_low_link(),other_vertex->get_index()));
        }
        sparse_edges_index++;
        
    }
    
    
    if(curr_vertex->get_low_link()==curr_vertex->get_index()){
        vertex* hold = st.top();
        st.pop();
        int hold_vertex_name;
        while( *hold != *curr_vertex && !st.empty()){
            hold_vertex_name = hold->get_vertex_name();
            vertices[hold_vertex_name]->set_low_link(curr_vertex->get_low_link());
            hold = st.top();
            st.pop();
        }
    }
}


/*
 * function that call's strong_connected_component function to find the components of the graph
 * @pram: pointer to array of pointers of sparse_matrix_element structs
 * @pram: the number of sprase_matrix_element's in sparse_graph
 * @pram: the number of vertices in the graph
 * @pram: stack used to perform tarjan's algorithm
 */

template <typename T>
std::vector<vertex*> graph_con_com(DenseMatrix<T>& sm){
    
    std::stack<vertex*> st;
    int num_vertices = sm.getNumberOfRows();
    std::vector<SparseElement<T> > sparse_form = sm.getSparseForm();
    std::vector<vertex*> vertices(num_vertices);
    int index=0;
    
    
    
    for(int j=0;j<num_vertices;j++)
    {
        vertices[j] = new vertex(j,-1);
    }
    
    for(int i=0;i < num_vertices;i++)
    {
           // compIdx[i]=-1;
    }
    
    for(int i=0;i < num_vertices;i++)
    {
        if( vertices[i]->get_low_link()==-1)
        {
            //printf("unset vertex %d has low link of %d\n", i, vertices[i].get_low_link());
            strong_com(sparse_form,num_vertices,&index,i,vertices,st);
        }
        else
        {
            //printf(" set vertex %d has low link of %d\n", i, vertices[i].get_low_link());
        }
    }
    
//    for(int j=0;j<num_vertices;j++)
//    {
//        std::cout<<"low link of: " <<  j << ": " << vertices[j].get_low_link() << std::endl;
//    }
    return vertices;
}


#endif
