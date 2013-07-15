//
//  Point.h
//  Graph Matching
//
//  Created by Ali Hajimirza on 7/10/13.
//  Copyright (c) 2013 Ali Hajimirza. All rights reserved.
//

#ifndef __Point__
#define __Point__

#include <iostream>
#include "MatrixExceptions.h"

class Point
{
private:
    void _copy(const Point&);
    
protected:
    int _x;
    int _y;
public:
    /**************
     *Constructors*
     **************/
    Point();
    Point(int, int);
    Point(const Point&);
    
    /************
     *Destructor*
     ************/
    virtual ~Point();
    
    /***********
     *ACCESSORS*
     ***********/
    int getX() const;
    int getY() const;

    /**********
     *MUTATORS*
     **********/
    void setX(int);
    void setY(int);
    
    /**********
     *OPERATORS*
     **********/
    bool operator>(const Point&) const;
    bool operator<(const Point&) const;
    bool operator==(const Point&) const;
    void operator=(const Point&);
    friend struct std::hash<Point>;
    friend std::ostream& operator<<(std::ostream&, const Point&);
};

template<>
struct std::hash<Point>
{
    std::size_t operator()(const Point& p) const
    {
        std::hash<int> int_hash;
        return ((51 + int_hash(p.getX())) * 51 + int_hash(p.getY()));
    }
};


#endif /* defined(__Graph_Matching__Point__) */
