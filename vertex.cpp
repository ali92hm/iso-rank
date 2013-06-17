#include "vertex.h"

/*vertex class default constructor
 */
vertex:: vertex(){
 this->vertex_name=-1;
 this->index=-1;
 this->low_link=100000;

}

/*vertex class constructor
 *parameter a is vertex_name b is the index
 */
vertex::vertex(int a,int b){
  this->vertex_name=a;
  this->index=b;
  this->low_link=1000000;
}

/*vertex class function
 *operator == that checks if this vertex is equal to parameter
 */
int vertex::operator==(const vertex &rhs){
  return this->vertex_name==rhs.vertex_name;
}

/*vertex class function
 *operator != that checks if this vertex is not equal to parameter
 */
int vertex::operator!=(const vertex &rhs){
  return 1-(this->vertex_name==rhs.vertex_name);

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
  this->low_link=b;
}

/*vertex class function
 *returns this vertex's index
 */
long vertex:: get_index(){
  return index;
}

/*vertex class function
 *returns this vertex's low link
 */
long vertex:: get_low_link(){
  return low_link;
}

/*vertex class function
 *returns this vertex's name
 */
int vertex:: get_vertex_name(){
return vertex_name;
}
