/*! \file contour.hpp
  \brief Finds contour points on an unstructured mesh
  \author Almog Yalinewich
 */

#ifndef CONTOUR_HPP
#define CONTOUR_HPP 1

#include <memory>
#include "main_loop_2d.hpp"
#include "index2filename.hpp"
#include "trigger.hpp"

class LocalContourCriterion
{
public:

  virtual std::pair<bool,Vector2D> operator()
  (const Edge& edge, const hdsim& sim) const = 0;

  virtual ~LocalContourCriterion(void);
};

class SequentialContour: public DiagnosticFunction
{
public:

  SequentialContour(Trigger* p_trigger,
		    Index2FileName* p_i2f,
		    LocalContourCriterion* p_lcc);

  void operator()(const hdsim& sim);

private:
  std::auto_ptr<Trigger> p_trigger_;
  int count_;
  std::auto_ptr<Index2FileName> p_i2f_;
  std::auto_ptr<LocalContourCriterion> p_lcc_;
};

#endif // CONTOUR_HPP
