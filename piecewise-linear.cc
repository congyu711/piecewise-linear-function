#include "piecewise-linear.hh"
#include <algorithm>
#include <iostream>
using namespace std;
const double eps=1e-9;
using plf::piecewise_linear_func;
using plf::line_segment;
double __doubleMax=std::numeric_limits<double>::max();
double plf::_intersection(const line_segment& f, const line_segment& g)
{
  if(f.a==g.a)  return __doubleMax;
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
piecewise_linear_func piecewise_linear_func::operator()(const piecewise_linear_func &g)
{
  piecewise_linear_func res;
  res.l=g.l; res.r=g.r;
  for(int i=0;i<g.lines.size();i++)
  {
    auto lv=g.lines[i].getval(g.lines[i].l), rv=g.lines[i].getval(g.lines[i].r);  // range
    auto fl=lower_bound(lines.begin(),lines.end(),std::min(lv,rv),[&](line_segment interval,double x){return interval.r<=x;})-lines.begin();
    auto fr=lower_bound(lines.begin(),lines.end(),std::max(lv,rv),[&](line_segment interval,double x){return interval.r<x;})-lines.begin();
    fr=std::min(lines.size()-1,(size_t)fr);  // if the range is too large.
    double _it;
    auto &line_g=g.lines[i];
    if(lv<rv)
    {
      _it=lv;
      for(int j=fl;j<=fr;j++)
      {
        auto &line_f=lines[j];
        if(line_f.r==_it)  continue;
        line_segment tmp((_it-line_g.b)/line_g.a,(std::min(line_f.r,rv)-line_g.b)/line_g.a,line_f.a*line_g.a,line_f.a*line_g.b+line_f.b);
        res.lines.push_back(tmp);
        _it=line_f.r;
      }
    }
    else if(lv>rv)
    {
      _it=lv;
      for(int j=fr;j>=fl;j--)
      {
        auto &line_f=lines[j];
        if(line_f.l==_it)  continue;
        line_segment tmp((_it-line_g.b)/line_g.a,(std::max(line_f.l,rv)-line_g.b)/line_g.a,line_f.a*line_g.a,line_f.a*line_g.b+line_f.b);
        res.lines.push_back(tmp);
        _it=line_f.l;
      }
    }
    else  // flat line
    {
      res.lines.push_back(line_segment(line_g.l,line_g.r,0,lines[fl].getval(lv)));
    }
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
    if(cross==__doubleMax||cross>min(_f.r,_g.r)||cross<left)
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

double piecewise_linear_func::operator()(double x)
{
  if(x>r||x<l) return __doubleMax;
  auto it=lower_bound(lines.begin(),lines.end(),x,[&](line_segment interval,double x){return interval.r<=x;});
  return (*it).getval(x);
}

piecewise_linear_func plf::_infconv(const piecewise_linear_func& f, const piecewise_linear_func& g)
{
  piecewise_linear_func res;
  res.l=std::min(f.l,g.l);
  // using slope_length=pair<double,double>;
  // vector<slope_length> slopes;
  double b0=0;
  // slopes vector should be sorted if f and g are convex functions.
  int i=0,j=0;
  bool __f=0,__g=0; // aux variable for overlapped segments
  // use one function only when another is not defined for some x;
  if(f.l<g.l)
  {
    while(f.lines[i].r<=g.l)
    {
      auto &c=f.lines[i];
      // slopes.push_back(make_pair(c.a,c.r-c.l));
      res.lines.push_back(c);
      i++;
    }
    // deal with the overlapped line segments
    if(g.l>=f.lines[i].l&& g.r<=f.lines[i].r)
    {
      // slopes.push_back(make_pair(f.lines[i].a, f.lines[i].r-g.l));
      res.lines.push_back(line_segment(f.lines[i].l,g.l,f.lines[i].a,f.lines[i].b));
      __f=1;
    }
  }
  else
  {
    while(g.lines[j].r<=f.l)
    {
      auto &c=g.lines[j];
      // slopes.push_back(make_pair(c.a,c.r-c.l));
      res.lines.push_back(c);
      j++;
    }
    // deal with the overlapped line segments
    if(f.l>=g.lines[j].l&& f.r<=g.lines[j].r)
    {
      // slopes.push_back(make_pair(g.lines[j].a, g.lines[j].r-f.l));
      res.lines.push_back(line_segment(g.lines[j].l,f.l,g.lines[j].a,g.lines[j].b));
      __g=1;
    }
  }

  for(;i<f.lines.size()||j<g.lines.size();)
  {
    // f
    if(i!=f.lines.size()&&(j==g.lines.size()||f.lines[i].a<=g.lines[j].a))
    {
      auto &c=f.lines[i];
      double len;
      if(__f) len=c.r-g.l, __f=0;
      else len=c.r-c.l;
      if(c.r-c.l!=0)  // slopes.push_back(make_pair(c.a,len));
      {
        auto& last=res.lines.back();
        res.lines.push_back(line_segment(last.r,last.r+len,c.a,last.getval(last.r)-c.a*last.r));
      }
      i++;
    }
    else if(j!=g.lines.size()&&(i==f.lines.size()||f.lines[i].a>=g.lines[j].a))
    {
      auto &c=g.lines[j];
      double len;
      if(__g) len=c.r-f.l, __g=0;
      else len=c.r-c.l;
      if(c.r-c.l!=0)  // slopes.push_back(make_pair(c.a,c.r-c.l));
      {
        auto& last=res.lines.back();
        res.lines.push_back(line_segment(last.r,last.r+len,c.a,last.getval(last.r)-c.a*last.r));
      }
      j++;
    }
  }
  res.r=res.lines.back().r;
  
  return res;
}