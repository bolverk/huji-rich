/*! \file SourceTerm3D.hpp
\brief Abstract class for source terms
\author Elad Steinberg
*/

#ifndef SOURCETERM3D_HPP
#define SOURCETERM3D_HPP 1
#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "../../misc/utils.hpp"
#include "computational_cell.hpp"
#include "conserved_3d.hpp"

//! \brief Abstract class for external forces
class SourceTerm3D
{
public:
	/*!
	\brief Calcualtes the change in conserved variables done on a cell from a source term
	\param tess The tessellation
	\param cells The hydrodynmic variables of the cell
	\param fluxes The hydrodynamic fluxes
	\param point_velocities Velocities of the mesh generating points
	\param t Time
	\param dt The time step
	\param tracerstickernames The names of the tracers and stickers
	\param extensives The updates extenesives, given as input and output.
	*/
	virtual void operator()(const Tessellation3D& tess,const vector<ComputationalCell3D>& cells,
		const vector<Conserved3D>& fluxes,const vector<Vector3D>& point_velocities,const double t,double dt,
		TracerStickerNames const& tracerstickernames,vector<Conserved3D> &extensives) const = 0;

	virtual double SuggestInverseTimeStep(void)const;

	virtual ~SourceTerm3D(void);
};

class ZeroForce3D : public SourceTerm3D
{
public:
	void operator()(const Tessellation3D& /*tess*/, const vector<ComputationalCell3D>& /*cells*/,
		const vector<Conserved3D>& /*fluxes*/, const vector<Vector3D>& /*point_velocities*/, const double /*t*/, 
		double /*dt*/,TracerStickerNames const& /*tracerstickernames*/, vector<Conserved3D> &/*extensives*/) const {}
};

#endif //SOURCETERM3D_HPP