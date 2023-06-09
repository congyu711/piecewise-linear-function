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
    double getval(double x) const
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
    
    // addition
    piecewise_linear_func operator+(const piecewise_linear_func& g);

    // composition f(g(x))
    // worst case running time is O(nm log n) where n and m are the numbers of breakpoints of f and g.
    piecewise_linear_func operator()(const piecewise_linear_func& g);

    // Find f(x). Assume the function is continuous.
    double operator()(double x);
  };
  // max
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

  // infimal convolution
  // f and g should be coutinuous piecewise linear convex functions.
  piecewise_linear_func _infconv(const piecewise_linear_func&, const piecewise_linear_func&);
}

#endif