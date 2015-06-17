/*! \file AlphaDisc.hpp
\brief Alpha model for disc viscosity, assumes center of disc is at origin
\author Elad Steinberg
*/

#ifndef ALPHADISC_HPP
#define ALPHADISC_HPP 1

#include "Viscosity.hpp"

//! \brief Gravitational acceleration due to a pointlike mass
class AlphaDisc : public Viscosity
{
public:
	/*!
	\brief Class constructor
	\param alpha The alpha parameter
	\param height_ratio The disk H/R ratio
	\param grad Reference to the spatial gradients
	*/
	AlphaDisc(double alpha, double height_ratio, SpatialReconstruction &grad);

private:
	
	double GetNu(Tessellation const& tess, PhysicalGeometry const& pg, vector<Primitive> const& cells, int point)const;

	double alpha_;
	double height_ratio_;
	SpatialReconstruction &grads_;
};

#endif // ALPHADISC_HPP
