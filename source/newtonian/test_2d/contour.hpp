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

//! \brief Recipe for creating a contour from the simulation
class LocalContourCriterion
{
public:

  /*! \brief Calculates the intersection of the contour with the line between neighboring mesh generating points
    \param edge Edge between neighbors
    \param sim Simulation
    \return A pair. First member is true if the contour intersect the line, and second member is the position of the intersection
   */
  virtual std::pair<bool,Vector2D> operator()
  (const Edge& edge, const hdsim& sim) const = 0;

  //! \brief Class destructor
  virtual ~LocalContourCriterion(void);
};

//! \brief Write contour files at consecutive times
class SequentialContour: public DiagnosticFunction
{
public:

  /*! \brief Class constructor
    \param p_trigger Trigger function
    \param p_i2f File name generator
    \param p_lcc Contour function
   */
  SequentialContour(Trigger* p_trigger,
		    Index2FileName* p_i2f,
		    LocalContourCriterion* p_lcc);

  void operator()(const hdsim& sim);

private:
  std::unique_ptr<Trigger> p_trigger_;
  int count_;
  std::unique_ptr<Index2FileName> p_i2f_;
  std::unique_ptr<LocalContourCriterion> p_lcc_;
};

#endif // CONTOUR_HPP
