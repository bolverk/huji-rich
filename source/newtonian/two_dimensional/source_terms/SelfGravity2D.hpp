/*! \file SelfGravity2D.hpp
  \brief Self gravity force module, just a simple Nbody tree code
\author Elad Steinberg
*/

#ifndef SELFGRAVITY_HPP
#define SELFGRAVITY_HPP 1

#include "../source_terms/ConservativeForce.hpp"
#include "../../../treecode/ANN.h"

//! \brief Self gravity force module, just a simple Nbody tree code
class SelfGravity: public Acceleration
{
public:
	/*!
	\brief Constructor
	\param OpenAngle The node openning criteria
	\param hFactor Smoothing length, if r<h F=1/(r^3+h^3) where h=hFactor*CellRadius
	*/
	SelfGravity(double OpenAngle=0.45,double hFactor=0.3);

	~SelfGravity(void);

  /*! \brief Calculates the acceleration due to self gravity
    \param tess Tessellation
    \param cell Fluid elements
    \param point Index of the mesh generating point
    \return Acceleration
    \todo Determine if this function is necessary
   */
  Vector2D Calculate(Tessellation const& tess,
		     vector<Primitive> const& cell,int point);

	Vector2D Calculate
		(Tessellation const& tess,
		vector<Primitive> const& cells,
		int point,vector<Conserved> const& fluxes,
		vector<Vector2D> const& point_velocity,
		HydroBoundaryConditions const& hbc,
		double t,
		double dt);

private:
	int _counter;
	int _length;
	double _OpenAngle;
	ANNkd_tree *_tree;
	ANNpointArray _treePoints;
	vector<double> masses;
	vector<double> h;
	double _Factor;
	SelfGravity(SelfGravity const& origin);
	SelfGravity& operator=(SelfGravity const& origin);
};

#endif // SELFGRAVITY_HPP
