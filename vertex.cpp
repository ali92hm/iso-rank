#include "vertex.h"

/*vertex class default constructor
 */
vertex:: vertex(){
 this->vertex_name=-1;
 this->index=-1;
 this->low_link=-1;
this->tarjan_flag=0;
}

/*vertex class constructor
 *parameter a is vertex_name b is the index
 */
vertex::vertex(int a,int b){
  this->vertex_name=a;
  this->index=b;
  this->low_link=-1;
    this->tarjan_flag=0;
}

vertex::~vertex()
{
    
}

/*vertex class function
 *operator == that checks if this vertex is equal to parameter
 */
bool vertex::operator==(const vertex &rhs){
  return this->vertex_name==rhs.vertex_name;
}

/*vertex class function
 *operator != that checks if this vertex is not equal to parameter
 */
bool vertex::operator!=(const vertex &rhs){
  return !(this->vertex_name==rhs.vertex_name);

}

/*vertex class function
 *operator = sets this vertex to parameter
 */
vertex* vertex::operator=(const vertex *rhs){
  this->vertex_name=rhs->vertex_name;
  this->index=rhs->index;
  this->low_link=rhs->low_link;
  return this;
}

/*vertex class function
 *operator = sets this vertex to parameter
 */
vertex& vertex::operator=(const vertex &rhs){
  this->vertex_name=rhs.vertex_name;
  this->index=rhs.index;
  this->low_link=rhs.low_link;
  return *this;
}

/*vertex class function
 *sets this vertex's index to a
 */
void vertex::set_index(long a){
    this->index=a;
  }

/*vertex class function
 *sets this vertex's low_link to b
 */
void vertex::set_low_link(long b){
    //std::cout<< "setting low link of "<<this->get_vertex_name()<<  "to:" << b << std::endl;
  this->tarjan_flag=1;
    this->low_link=b;
}

/*vertex class function
 *returns this vertex's index
 */
long vertex:: get_index(){
  return this->index;
}

/*vertex class function
 *returns this vertex's low link
 */
long vertex:: get_low_link(){
  return this->low_link;
}

/*vertex class function
 *returns this vertex's name
 */
int vertex:: get_vertex_name(){
return this->vertex_name;
}

int vertex:: get_tarjan_flag(){
    return this->tarjan_flag;
}

/*
 * takes an array of vertex objects vertices and an integer
 * signifying a component of the graph. The function returns
 * an integer array where vertices that don't belong to the
 * component are masked out.
 */
int* component_mask(vertex* vertices, int component, int number_nodes){
    int* comp_mask= new int[number_nodes];//(int *)malloc(sizeof(int) * number_nodes);
    
    for(int i=0;i<number_nodes;i++){
        comp_mask[i]=0;
    }
    
    for(int i=0;i<number_nodes;i++){
        vertex* curr_vertex= &vertices[i];
        if((*curr_vertex).get_low_link()==component){
            comp_mask[i]= 1;
        }
    }
    
    return comp_mask;
}
