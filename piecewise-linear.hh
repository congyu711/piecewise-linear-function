#ifndef _PIECEWISE_LINEAR_
#define _PIECEWISE_LINEAR_

#include <vector>
#include <limits>
// #include <algorithm>

namespace plf // piecewise linear function
{

  class line_segment
  {
  public:
    double l, r;
    double a, b; // y=ax+b, x \in [l,r]
    line_segment(double _l,double _r,double _a,double _b):l(_l),r(_r),a(_a),b(_b){}
    line_segment():l(0),r(0),a(0),b(0){}
    double getval(double x)
    {
      return a*x+b;
    }
  };
  double _intersection(const line_segment& f, const line_segment& g);

  class piecewise_linear_func
  {
  public:
    double l,r;
    std::vector<line_segment> lines;   // assume line segments in this vector are always ordered and disjoint
                                  // y=0 for undefined intervals.
    piecewise_linear_func operator+(const piecewise_linear_func& g);
  };
  piecewise_linear_func _max(const piecewise_linear_func& f, const piecewise_linear_func& g);

  template<typename T>
  T max(const T& first)
  {
    return first;
  }
  template<typename T, typename... Args>
  T max(const T& first, Args... args)
  {
    return _max(first,max(args...));
  }

}

#endif