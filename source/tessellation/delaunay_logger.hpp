/*! \file delaunay_logger.hpp
  \brief A debugging function for the Delaunay triangulation
  \author Almog Yalinewich
 */

#ifndef DELAUNAY_OUTPUT_HPP
#define DELAUNAY_OUTPUT_HPP 1

#include <vector>
#include "geometry.hpp"
#include "facet.hpp"

//! \brief Diagnostic functions for Delaunay triangulation
namespace delaunay_loggers {
  using namespace std;

  //! \brief Diagnostic class for Delaunay triangulation
  class DelaunayLogger
  {
  public:

    //! \brief Class constructor
    DelaunayLogger(void);

    //! \brief Dumps output
    virtual void output(vector<Vector2D> const& cor,
			vector<facet> const& f);

    virtual ~DelaunayLogger(void);
  };

  //! \brief Writes data to a binary file
  class BinaryLogger: public DelaunayLogger
  {
  public:

    /*! \brief Class constructor
      \param file_name Name of output file
     */
    explicit BinaryLogger(string const& file_name);

    void output(vector<Vector2D> const& cor,
		vector<facet> const& f);

  private:
    const string file_name_;
  };
}

#endif // DELAUNAY_OUTPUT_HPP
