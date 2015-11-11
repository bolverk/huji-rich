/*! \file voronoi_logger.hpp
  \brief A debugging function for the Voronoi tessellation
  \author Almog Yalinewich
 */

#ifndef VORONOI_LOGGER_HPP
#define VORONOI_LOGGER_HPP 1

#include <string>
#include <vector>
#include "Edge.hpp"

class VoronoiMesh;
class Tessellation;

//! \brief Diagnostic functions for Voronoi tessellation
namespace voronoi_loggers
{
  using namespace std;

  //! \brief Abstract class for voronoi tessellation
  class VoronoiLogger
  {
  public:

    //! \brief Class constructor
    VoronoiLogger(void);

    /*! \brief Outputs information from the Voronoi tessellation
      \param v Refernce to Voronoi tessellation
     */
    virtual void output(const VoronoiMesh& v);

	 /*! \brief Outputs information from a tessellation
      \param v Refernce to the tessellation
     */
    virtual void output(const Tessellation& v);

    virtual ~VoronoiLogger(void);
  };

  //! \brief Writes data to a binary file
  class BinLogger: public VoronoiLogger
  {
  public:

    /*! \brief Class constructor
      \param file_name Name of output file
     */
    explicit BinLogger(const std::string& file_name);

    void output(const VoronoiMesh& v);

	void output(const Tessellation& v);

	/*! \brief Reads the output information from the Voronoi tessellation
      \param location Name of output file
	  \return The mesh points
     */
	vector<Vector2D> read(string location);

  private:
    const std::string file_name_;
  };
}

#endif // VORONOI_LOGGER_HPP
