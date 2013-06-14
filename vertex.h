


#ifndef _vertex_h
#define _vertex_h


class vertex{
  int vertex_name,index,low_link;
public:
  vertex();
  vertex (int, int);
  vertex* operator=(const vertex *rhs);
  vertex& operator=(const vertex &rhs);
  int operator!=(const vertex &rhs);
  int operator==(const vertex &rhs);
  void set_index(int);
  void set_low_link(int);
  int get_index();
  int get_low_link();
  int get_vertex_name();
};

#endif