


#ifndef _vertex_h
#define _vertex_h

#include <iostream>

class vertex{
    int vertex_name;
    long index;
    long low_link;
    int tarjan_flag;
public:
  vertex();
  vertex (int, int);
  vertex* operator=(const vertex *rhs);
  vertex& operator=(const vertex &rhs);
  int operator!=(const vertex &rhs);
  int operator==(const vertex &rhs);
  void set_index(long);
  void set_low_link(long);
  long get_index();
  long get_low_link();
  int get_vertex_name();
  int get_tarjan_flag();
};

#endif