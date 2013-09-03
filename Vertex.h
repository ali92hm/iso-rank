/***********************************************************************************
 * Vertex Class. Used in Tarjan's Algorithm to keep track of which vertex belongs  *
 * to which component.                                                             *
 ***********************************************************************************/


#ifndef _Vertex_h
#define _Vertex_h

class vertex
{
    int vertex_name;
    long index;
    long low_link;
    int tarjan_flag;
    
public:
    vertex();
    vertex (int, int);
    virtual ~vertex();
    vertex* operator=(const vertex *rhs);
    vertex& operator=(const vertex &rhs);
    bool operator!=(const vertex &rhs);
    bool operator==(const vertex &rhs);
    void set_index(long);
    void set_low_link(long);
    long get_index();
    long get_low_link();
    int get_vertex_name();
    int get_tarjan_flag();
};
#endif