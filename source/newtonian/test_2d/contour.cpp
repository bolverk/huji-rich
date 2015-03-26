#include <boost/foreach.hpp>
#include "contour.hpp"
#include <fstream>
#include "../../misc/simple_io.hpp"
#include "../../misc/hdf5_utils.hpp"
#include "../../misc/lazy_list.hpp"

LocalContourCriterion::~LocalContourCriterion(void) {}

namespace 
{
  vector<Vector2D> calc_contour_points(const hdsim& sim,
				       const LocalContourCriterion& lcc)
  {
    vector<Vector2D> res;
    const Tessellation& tess = sim.getTessellation();
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
(Trigger* p_trigger,
 Index2FileName* p_i2f,
 LocalContourCriterion* p_lcc):
  p_trigger_(p_trigger),
  count_(0),
  p_i2f_(p_i2f),
  p_lcc_(p_lcc) {}

namespace {

  class ComponentExtractor: public LazyList<double>
  {
  public:

    ComponentExtractor(const vector<Vector2D>& source,
		       double Vector2D::* component):
      source_(source), component_(component) {}

    size_t size(void) const
    {
      return source_.size();
    }

    double operator[](size_t i) const
    {
      return source_[i].*component_;
    }

  private:
    const vector<Vector2D>& source_;
    double Vector2D::* component_;
  };

  void hdf5_write(const string& output_file,
		  const vector<Vector2D>& data,
		  const double time)
  {
    if(!data.empty())
      (HDF5Shortcut(output_file))
	("time",vector<double>(1,time))
	("x",serial_generate(ComponentExtractor(data, &Vector2D::x)))
	("y",serial_generate(ComponentExtractor(data, &Vector2D::y)));
  }
}

void SequentialContour::operator()(const hdsim& sim)
{
  if((*p_trigger_.get())(sim)){
    hdf5_write((*p_i2f_.get())(count_),
	       calc_contour_points(sim,*p_lcc_.get()),
	       sim.getTime());
    ++count_;
  }
}
