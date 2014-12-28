#include <boost/foreach.hpp>
#include "contour.hpp"
#include <fstream>
#include "../../misc/simple_io.hpp"

namespace 
{
  vector<Vector2D> calc_contour_points(const hdsim& sim,
				       const LocalContourCriterion& lcc)
  {
    vector<Vector2D> res;
    const Tessellation& tess = sim.GetTessellation();
    BOOST_FOREACH(const Edge& edge,tess.getAllEdges())
      {
	const std::pair<int,Vector2D> temp = lcc(edge,sim);
	if(temp.first)
	  res.push_back(temp.second);
      }
    return res;
  }
}

SequentialContour::SequentialContour
(std::auto_ptr<Trigger> p_trigger,
 std::auto_ptr<Index2FileName> p_i2f,
 std::auto_ptr<LocalContourCriterion> p_lcc):
  p_trigger_(p_trigger),
  count_(0),
  p_i2f_(p_i2f),
  p_lcc_(p_lcc) {}

namespace {
  void ascii_write(const string& output_file,
		   const vector<Vector2D>& data)
  {
    std::ofstream f(output_file.c_str());
    BOOST_FOREACH(const Vector2D& datum, data)
      {
	f << datum.x << " " << datum.y << "\n";
      }
    f.close();
  }
}

void SequentialContour::operator()(const hdsim& sim)
{
  if((*p_trigger_.get())(sim)){
    ascii_write((*p_i2f_.get())(count_),calc_contour_points(sim,*p_lcc_.get()));
    write_number(sim.GetTime(),"timestamp_"+(*p_i2f_.get())(count_));
    ++count_;
  }
}
