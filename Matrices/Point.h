/*********************************************************************************
 * This class represent a point with a x and y coordinate.                       *
 *********************************************************************************/

#ifndef __Point__
#define __Point__

#include <iostream>
#include "MatrixExceptions.h"

/*
 * Point class definition and method declaration.
 */
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
    friend std::ostream& operator<<(std::ostream&, const Point&);
};
#endif
