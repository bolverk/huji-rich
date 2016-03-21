#include <boost/foreach.hpp>
#include "duffell.hpp"

Duffell::Duffell
(const double alpha,
 const int iter):
  alpha_(alpha),
  iter_(iter) {}

namespace {
  vector<Vector2D> extract_velocity
  (const vector<ComputationalCell>& cells,size_t N)
  {
    vector<Vector2D> res(N);
    for(size_t i=0;i<res.size();++i)
      res.at(i) = cells.at(i).velocity;
    return res;
  }

  Vector2D neighbor_average
  (const Tessellation& tess,
   const vector<Vector2D>& previous,
   const size_t i)
  {
    Vector2D res(0,0);
    double circumference = 0;
    const vector<int> edge_indices = 
      tess.GetCellEdges(static_cast<int>(i));
    BOOST_FOREACH(int index, edge_indices){
      const Edge edge = tess.GetEdge(index);
      const int other = 
	edge.neighbors.first + 
	edge.neighbors.second -
	static_cast<int>(i);
      const double edge_length = 
	abs(edge.vertices.first-
	    edge.vertices.second);
      circumference += edge_length;
      if(other<tess.GetPointNo())
	res += edge_length*
	  previous.at(static_cast<size_t>(other));
    }
    return (1.0/circumference)*res;
  }

  vector<Vector2D> smooth
  (const Tessellation& tess,
   const vector<Vector2D>& previous,
   const double alpha)
  {
    vector<Vector2D> res(previous.size());
    for(size_t i=0;i<res.size();++i)
      res.at(i) = previous.at(i)*(1-alpha)+
	alpha*neighbor_average(tess,previous,i);
    return res;
  }
}

vector<Vector2D> Duffell::operator()
(const Tessellation& tess,
 const vector<ComputationalCell>& cells,
 double /*time*/,TracerStickerNames const& /*tracerstickernames*/) const
{
  vector<Vector2D> res = extract_velocity(cells,static_cast<size_t>(tess.GetPointNo()));
  for(int i=0;i<iter_;++i)
    res = smooth(tess,res,alpha_);
  return res;
}
