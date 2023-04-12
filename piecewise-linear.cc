#include "piecewise-linear.hh"
#include <algorithm>
using namespace plf;
using namespace std;
const double eps=1e-9;

piecewise_linear_func piecewise_linear_func::operator+(const piecewise_linear_func &g)
{
  piecewise_linear_func res;
  res.l=std::min(l,g.l);
  res.r=std::max(r,g.r);
  vector<double> bkpts,g_bkpts; // breakpoints
  vector<vector<line_segment>::const_iterator> ptrs,g_ptrs; // pointers
  bkpts.push_back(res.l);g_bkpts.push_back(res.l);
  for(auto it=lines.cbegin();it!=lines.cend();it++)
  {
    auto e=*it;
    if(e.l!=bkpts.back()) // there is an undefined interval
    {
      bkpts.push_back(e.l);
      ptrs.push_back(lines.cend());
    }
    bkpts.push_back(e.r);
    ptrs.push_back(it);
  }
  if(bkpts.back()!=res.r)
  {
    bkpts.push_back(res.r);
    ptrs.push_back(lines.cend());
  }
  for(auto it=g.lines.cbegin();it!=g.lines.cend();it++)
  {
    auto e=*it;
    if(e.l!=g_bkpts.back()) // there is an undefined interval
    {
      g_bkpts.push_back(e.l);
      g_ptrs.push_back(g.lines.cend());
    }
    g_bkpts.push_back(e.r);
    g_ptrs.push_back(it);
  }
  if(g_bkpts.back()!=res.r)
  {
    g_bkpts.push_back(res.r);
    g_ptrs.push_back(g.lines.cend());
  }
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