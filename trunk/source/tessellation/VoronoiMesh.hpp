/*! \file VoronoiMesh.hpp
  \brief Voronoi tessellation
  \author Elad Steinberg
 */

#ifndef VORONOIMESH_HPP
#define VORONOIMESH_HPP 1

#define _USE_MATH_DEFINES
#include "Delaunay.hpp"
#include "tessellation.hpp"
#include "../misc/utils.hpp"
#include "../newtonian/two_dimensional/HydroBoundaryConditions.hpp"
#include "../newtonian/two_dimensional/RefineStrategy.hpp"

//! \brief Voronoi tessellation class
class VoronoiMesh : public Tessellation
{
public:
	/*!
	\brief Removes cells from a Voronoi tessellation
	\param V The Voronoi to change
	\param ToRemove The indeces of the cells to remove
	\param VolIndex The indeces of the neighbors of each removed cell, the input of the vector is overwritten with the output
	\param Volratio How much volume relative to the removed cell's volume did each of the neighbors recieve, the input of the vector is overwritten with the output
	*/
	friend void Remove_Cells(VoronoiMesh &V,vector<int> &ToRemove,
		vector<vector<int> > &VolIndex,vector<vector<double> > &Volratio);
	/*!
	\brief Refines cells
	\param V The Voronoi tesselllation
	\param ToRefine The list of points to split and refine
	\param alpha The relative seperation between the splitted points (in units of the cell radius)
	\param directions The directions to give the splitted points (can be given empty)
	\param PeriodicBoundary Flag whether the grid is periodic or not
	*/
	friend void Refine_Cells(VoronoiMesh &V,vector<int> const& ToRefine,
		double alpha,vector<Vector2D> const& directions,bool PeriodicBoundary);

	vector<Vector2D> calc_edge_velocities(HydroBoundaryConditions const* hbc,
		vector<Vector2D> const& point_velocities,double time)const;

	Vector2D CalcFaceVelocity(Vector2D wl, Vector2D wr,Vector2D rL,
		Vector2D rR,Vector2D f)const;

	void RemoveCells(vector<int> &ToRemove,vector<vector<int> > &VolIndex,
		vector<vector<double> > &Volratio);

	void RefineCells(vector<int> const& ToRefine,
		vector<Vector2D> const& directions,double alpha=1e-3);

	Tessellation* clone(void)const;

	vector<int> GetNeighbors(int index)const;

	int GetOriginalIndex(int point) const;
	
	void Initialise(vector<Vector2D> const& points, OuterBoundary const* _bc);

	//! \brief Class default constructor.
	VoronoiMesh(void);

	/*!
		\brief Copy constructor
		\param other The Voronoi to copy
	*/
	VoronoiMesh(VoronoiMesh const& other);

	/*! 
		\brief Update voronoi tessellation
		\param points The new points
   */
	void Update(vector<Vector2D> const& points);

	~VoronoiMesh(void);
	/*! \brief Get Total number of mesh generating points
    \return Number of mesh generating points
   */
  int GetPointNo(void) const;

  /*! \brief Returns Position of mesh generating point
    \param index Mesh generating point index
    \return Position of mesh generating point
   */
  Vector2D GetMeshPoint(int index) const;

  /*! \brief Returns the total number of faces
    \return Total number of faces
   */
  int GetTotalSidesNumber(void) const;

  /*! \brief Returns edge (interface between cells)
    \param index Face index
    \return Interface between cells
   */
  Edge const& GetEdge(int index) const;

  /*! \brief Returns the effective width of a cell
    \param index Cell index
    \return Effective cell width
   */
  double GetWidth(int index) const;

  /*! \brief Returns the volume of a cell
    \param index Cell index
    \return Cell volume 
   */
  double GetVolume(int index) const;
  
//! \brief Returns a reference to a vector<int> containing the indexes of the edges related to a cell. \param index The index of the cell. \returns The reference to vector<int> 
  vector<int>const& GetCellEdges(int index) const;
  
/*! \brief Returns Position of Cell's CM
    \param index Mesh generating point index (the cell's index)
    \return Position of CM
   */
  Vector2D const& GetCellCM(int index) const;

  vector<Vector2D>& GetMeshPoints(void);

  /*! \brief Outputs the grid data
  \param filename The path to the output file
  */
  void output(string filename);
	
  bool NearBoundary(int index) const;

private:
	OuterBoundary const* bc;
	double eps;
	vector<Edge> edges;
	vector<vector< int> > added;
	vector<Vector2D> CM;
	vector<vector<int> > mesh_vertices; // Which edges does each mesh point have
	Delaunay *Tri;
	Vector2D get_center(int facet);
	bool legal_edge(Edge *e);
	void build_v(void);//Builds the voronoi mesh
	VoronoiMesh& operator=(const VoronoiMesh& origin);
	Vector2D CalcCellCM(int index) const;
};


#endif // VORONOIMESH_HPP
