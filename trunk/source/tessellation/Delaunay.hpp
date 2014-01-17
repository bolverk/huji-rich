/*! \file Delaunay.hpp
  \brief Delaunay triangulation
  \author Elad Steinberg
 */

#ifndef DELAUNAY_HPP
#define DELAUNAY_HPP 1
#include "facet.hpp"
#include <algorithm>
#include <vector>
#include <stack>
#include <fstream>
#include <climits>
#include "../../source/tessellation/geometry.hpp"
#include "tessellation.hpp"
#include "../newtonian/two_dimensional/OuterBoundary.hpp"
#include "../treecode/ANN.h"			
#include "../treecode/ANNperf.h"
#include "../misc/universal_error.hpp"
#include "HilbertOrder.hpp"
#include "geotests.hpp"
#include "../misc/utils.hpp"

using namespace std;
/*! \brief The Delaunay data structure. Gets a set of points and constructs the Delaunay tessellation.
\author Elad Steinberg
*/
class Delaunay
{
private:
	enum Sides{RIGHT,UP,LEFT,DOWN,LU,LD,RU,RD};
	ANNkd_tree *tree;
	ANNpointArray treePoints;
	vector<int> NewPointIndex;
	int lastFacet; //last facet to be checked in Walk
	bool CalcRadius;

	class DataOnlyForBuild
	{
	public:		
		DataOnlyForBuild();
		DataOnlyForBuild(DataOnlyForBuild const& other);
		DataOnlyForBuild& operator=(DataOnlyForBuild const& other);
		vector<int> insert_order;
		vector<vector<char> > copied;
		vector<int> BoundaryCandidates;
	};

	vector<double> radius;
	bool PointWasAdded;
	double BIG;
	int last_facet_added;
	vector<facet> f;	
	vector<Vector2D> cor;
	vector<Vector2D> totalCor;
	int length,olength;
	int location_pointer;
	int last_loc;
	OuterBoundary const* bc;

	void BuildTree(DataOnlyForBuild &data);
	bool BoundaryCandidate(double maxR,int point);
	double GetInitialSearchR(int index,double &maxradius);
	int FindNextOuterFacet(int facet,int point,int &nextpoint);
	void GetOuterPoints(vector<vector<int> > &OuterPoints);
	void GetOuterPoints2(vector<vector<int> > &OuterPoints);
	void DuplicatePoints(DataOnlyForBuild &data);
	void add_point(int index);
	void flip(int i,int j);
  //	int find_index(facet *fc,int i) const;
	void find_diff(facet *f1,facet *f2,int*) const;
	int Walk(int point);
	void CheckInput(void);
	int AddSphere(double minR,double MaxR,ANNpoint &queryPt,Sides side,
		bool partial,DataOnlyForBuild &data);
	bool GetQueryPointReflective(Vector2D &point,Sides side,double SearchR,int p);
	bool GetQueryPointPeriodic(Vector2D &point,Sides side,double SearchR,int p);
	void AddBoundaryPoint(int pointToAdd,Sides side,bool parital,DataOnlyForBuild
		 &data);
	double CalculateRadius(int facet);
	bool IsOuterFacet(int facet);
	int FindPointInFacet(int facet,int point);
	double FindMaxRadius(int point);
	void FindContainingTetras(int StartTetra,int point,vector<int> &tetras);
	void BuildBoundary(DataOnlyForBuild &data);
	double ReflectDistance(double SearchR,Sides side,int point,int numOfPoints,int &closest);

	//void DuplicatePointsLimited2(void);
	//void AddBigTirangle(int point);
	//vector<double> GetNonLocalRadii(int point,vector<Vector2D> &centers);
	//Vector2D CalCircleCenter(int facet);
	//void BuildBoundaryLimited(void);
	//void DuplicatePointsLimited(void);


  Delaunay& operator=(const Delaunay& origin);
public:
	/*! \brief Changes the cor olength
	\param n The new length;
	*/
	void ChangeOlength(int n);

	/*! \brief Changes the cor length
	\param n The new length
	*/
	void Changelength(int n);

	/*! \brief Allows to change the NewPointIndex
	\return Refrence to the NewPointIndex vector
	*/
	vector<int>& ChangeNewPointIndex(void);

	/*! \brief Allows to change the cor
	\return Refrence to the cor vector
	*/
	vector<Vector2D>& ChangeCor(void);

	//! \brief Class constructor.
	Delaunay(void);

	/*!
		\brief Copy constructor
		\param other The Triangulation to copy
	*/
	Delaunay(Delaunay const& other);

	//! \brief Default destructor.
	~Delaunay(void);
	
	/*! \brief Builds the Delaunay tessellation.  
	\param vp A refrence to a vector of points to be added. 
	\param _bc The location of the boundary.
	*/
	void build_delaunay(vector<Vector2D>const& vp,OuterBoundary const *_bc);
	
	//! \brief Dumps the Delaunay tessellation into a binary file.
	void output(void);

	/*! \brief Returns a facet.
	\param index The faet to return.
	\returns A pointer to the selected facet.
	*/
	facet* get_facet(int index);

	/*! \brief Returns a coordinate of a vertice. 
	\param Facet The index of the facet to check. 
	\param vertice The index of the vertice in the facet. 
	\param dim If dim=0 returns the x-coordinate else returns the y-coordinate. 
	\returns The chosen coordinate.
	*/
	double get_facet_coordinate(int Facet,int vertice, int dim);

	/*! \brief Returns a point.
	\param index The index of the point. 
	\returns The chosen point.
	*/
	Vector2D get_point(int index) const;

	/*! \brief Returns a coordinate.
	\param index The index of the point. 
	\param dim If dim=0 returns the x-coordinate else returns the y-coordinate. 
	\returns The chosen coordinate.
	*/
	double get_cor(int index,int dim);

	/*! \brief Returns the number of facets. 
	\returns The number of facets.
	*/
	int get_num_facet(void);

	/*! \brief Returns the number of points
	\returns The number of points.
	*/
	int get_length(void) const;

	/*! \brief Returns the last location, a number used to identify the fact that the neighbor of a facet is empty. 
	\returns The last location.
	*/
	int get_last_loc(void) const;

	/*! \brief Change Mesh point. 
	\param index The index of the point to change. 
	\param p The new point to set.
	*/
	void set_point(int index, Vector2D p);

	/*! \brief Returns the area of the triangle. Negative result means the triangle isn't right handed. 
	\param index The index to the facet 
	\return The area
	*/
	double triangle_area(int index);

	/*!
	\brief Updates the triangulation
	\param points The new set of points
	*/
	void update(const vector<Vector2D>& points);
	
	/*!
		\brief Returns the original index of the duplicated point in Periodic Boundary conditions
		\param NewPoint The index of the duplicated point
		\return The original index
	*/
	int GetOriginalIndex(int NewPoint) const;

	/*!
		\brief Returns the original length of the points (without duplicated points)
		\return The original length
	*/
	int GetOriginalLength(void) const;

	/*!
		\brief Returns a refrence to the points
		\return Refrence to the points
	*/
	vector<Vector2D>& GetMeshPoints(void);

	/*! \brief Returns the length of all the points (included duplicated)
	\return The length of all of the points
	*/
	int GetTotalLength(void);

	/*! \brief return the facet's radius
	\param facet The facet to check
	\return The facet's radius
	*/
  double GetFacetRadius(int facet) const;
};
#endif //DELAUNAY_HPP
