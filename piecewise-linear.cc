#include "piecewise-linear.hh"
#include <algorithm>

using namespace std;
const double eps=1e-9;
using plf::piecewise_linear_func;
using plf::line_segment;
double plf::_intersection(const line_segment& f, const line_segment& g)
{
  if(f.a==g.a)  return std::numeric_limits<double>::max();
  return (g.b-f.b)/(f.a-g.a);
}

void linesegment2breakpoints(const piecewise_linear_func &g, double l, double r, 
                             vector<double>& bkpts, vector<vector<line_segment>::const_iterator>& ptrs)
{
  bkpts.clear();ptrs.clear();
  bkpts.push_back(l);
  for(auto it=g.lines.cbegin();it!=g.lines.cend();it++)
  {
    auto e=*it;
    if(e.l!=bkpts.back()) // there is an undefined interval
    {
      bkpts.push_back(e.l);
      ptrs.push_back(g.lines.cend());
    }
    bkpts.push_back(e.r);
    ptrs.push_back(it);
  }
  if(bkpts.back()!=r)
  {
    bkpts.push_back(r);
    ptrs.push_back(g.lines.cend());
  }
}

piecewise_linear_func piecewise_linear_func::operator+(const piecewise_linear_func &g)
{
  piecewise_linear_func res;
  res.l=std::min(l,g.l);
  res.r=std::max(r,g.r);
  vector<double> bkpts,g_bkpts; // breakpoints
  vector<vector<line_segment>::const_iterator> ptrs,g_ptrs; // pointers
  linesegment2breakpoints(ref(*this),res.l,res.r,bkpts,ptrs);
  linesegment2breakpoints(g,res.l,res.r,g_bkpts,g_ptrs);

  double left=res.l;
  for(int i=1,j=1;i<bkpts.size()&&j<g_bkpts.size();)
  {
    line_segment tmp;
    tmp.l=left;
    if(ptrs[i-1]!=lines.cend()) tmp.a+=ptrs[i-1]->a,tmp.b+=ptrs[i-1]->b;
    if(g_ptrs[j-1]!=g.lines.cend()) tmp.a+=g_ptrs[j-1]->a,tmp.b+=g_ptrs[j-1]->b;
    if(bkpts[i]<g_bkpts[j])
    {
      tmp.r=bkpts[i];
      i++;
    }
    else if(bkpts[i]>g_bkpts[j])
    {
      tmp.r=g_bkpts[j];
      j++;
    }
    else  // bkpts[i]==g_bkpts[j]
    {
      tmp.r=bkpts[i];
      i++,j++;
    }
    left=tmp.r;
    res.lines.push_back(tmp);
  }

  return res;
}

piecewise_linear_func plf::_max(const piecewise_linear_func& f, const piecewise_linear_func& g)
{
  piecewise_linear_func res;
  res.l=std::min(f.l,g.l);
  res.r=std::max(f.r,g.r);
  vector<double> f_bkpts,g_bkpts; // breakpoints
  vector<vector<line_segment>::const_iterator> f_ptrs,g_ptrs; // pointers
  linesegment2breakpoints(f,res.l,res.r,f_bkpts,f_ptrs);
  linesegment2breakpoints(g,res.l,res.r,g_bkpts,g_ptrs);
  double left=res.l;
  for(int i=1,j=1;i<f_bkpts.size()&&j<g_bkpts.size();)
  {
    line_segment tmp,_f,_g;
    tmp.l=left;
    if(f_ptrs[i-1]!=f.lines.cend()) _f=*f_ptrs[i-1];
    else _f.r=f_bkpts[i];
    if(g_ptrs[j-1]!=g.lines.cend()) _g=*g_ptrs[j-1];
    else _g.r=g_bkpts[j];
    _f.l=_g.l=left;
    // compute the intersection
    double cross=_intersection(_f,_g);
    if(_f.getval(left)<_g.getval(left)) tmp.a=_g.a,tmp.b=_g.b;
    else tmp.a=_f.a,tmp.b=_f.b;
    // no valid intersection.
    if(cross==numeric_limits<double>::max()||cross>min(_f.r,_g.r)||cross<left)
    {
      tmp.r=min(f_bkpts[i],g_bkpts[j]);
      res.lines.push_back(tmp);
    }
    else
    {
      tmp.r=cross;
      res.lines.push_back(tmp);
      // add another line segment
      double __a=(tmp.a==_f.a?_g.a:_f.a);
      double __b=(tmp.a==_f.a?_g.b:_f.b);
      res.lines.push_back(line_segment(cross,min(f_bkpts[i],g_bkpts[j]),__a,__b));
    }

    left=min(f_bkpts[i],g_bkpts[j]);
    if(f_bkpts[i]<g_bkpts[j])
      i++;
    else if(f_bkpts[i]>g_bkpts[j])
      j++;
    else  // bkpts[i]==g_bkpts[j]
      i++,j++;

  }
  return res;
}