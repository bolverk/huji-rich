/*! \brief Base class for disc tessellation
\author Elad Steinberg
*/

#ifndef DISC_TESS
#define DISC_TESS 1

#include "../tessellation/tessellation.hpp"
#include "../misc/utils.hpp"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "../newtonian/two_dimensional/HydroBoundaryConditions.hpp"
#define PI 3.14159265359

class DiscTess: public Tessellation
{
private:
	double _MaxR,_MinR;
	vector<double> _r;
	vector<Vector2D> cor;
	vector<vector<int> > _mesh_vertices;
	vector<Edge> _edges;
	int PointNum;

public:
	/*!
	\brief Class constructor
	\param MaxR The disc maximum radius
	\param MinR The disc minimum radius
	\param R The radial distances of the annuli (sorted)
	*/
	DiscTess(double MaxR,double MinR,vector<double> R);

	~DiscTess(void);

	void Initialise(vector<Vector2D> const& points,OuterBoundary const* bc=0);

	void Update(vector<Vector2D> const& points);

	int GetPointNo(void) const;

	Vector2D GetMeshPoint(int index) const;

	Vector2D GetCellCM(int index) const;

	int GetTotalSidesNumber(void) const;

	Edge const& GetEdge(int index) const;

	double GetWidth(int index) const;

	double GetVolume(int index) const;

	vector<int>const& GetCellEdges(int index) const;

	int GetOriginalIndex(int point) const;

	vector<Vector2D>& GetMeshPoints(void);

	vector<int> GetNeighbors(int index)const;

	Tessellation* clone(void) const;

	bool NearBoundary(int index) const;

	void RemoveCells(vector<int> &ToRemovevector,vector<vector<int> > &VolIndex,
		vector<vector<double> > &Volratio);

	void RefineCells(vector<int> const& ToRefine,vector<Vector2D> const&
		directions,double alpha);

	/*! \brief Outputs the grid data
	\param filename The path to the output file
	*/
	void output(string filename);

	Vector2D CalcFaceVelocity(Vector2D wl, Vector2D wr,Vector2D rL,
		Vector2D rR,Vector2D f)const;

	vector<Vector2D> calc_edge_velocities(HydroBoundaryConditions const* hbc,
		vector<Vector2D> const& point_velocities,double time)const;

};

#endif // DISC_TESS