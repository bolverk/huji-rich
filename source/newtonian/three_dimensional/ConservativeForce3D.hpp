/*! \file ConservativeForce3D.hpp
\brief Abstract class for conservative force's acceleration
\author Elad Steinberg
*/

#ifndef CONSFORCE3D_HPP
#define CONSFORCE3D_HPP 1

#include "SourceTerm3D.hpp"

//! \brief Physical acceleration
class Acceleration3D
{
public:
	/*!
	\brief Calculates the acceleration that the cells feel
	\param tess The tessellation
	\param cells The primitive cells
	\param fluxes The vector of the fluxes
	\param time The simulation time
	\param tracerstickernames The names of the tracers and stickers
	\param acc The calculated acceleration, given as output
	*/
	virtual void operator()(const Tessellation3D& tess,const vector<ComputationalCell3D>& cells,
			const vector<Conserved3D>& fluxes,const double time,TracerStickerNames const& tracerstickernames,
		vector<Vector3D> &acc) const = 0;

	virtual ~Acceleration3D(void);
};

/*! \brief Class for conservative forces
\author Elad Steinberg
*/
class ConservativeForce3D : public SourceTerm3D
{
public:
	/*! \brief Class constructor
	\param acc The acceleration force
	\param mass_flux TO include the mass flux term in the calculation or not (eq.  94 or eq. 82 in Arepo)
	*/
	explicit ConservativeForce3D(const Acceleration3D& acc,bool mass_flux);

	/*!
	\brief Class destructor
	*/
	~ConservativeForce3D(void);

	void operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells,
		const vector<Conserved3D>& fluxes, const vector<Vector3D>& point_velocities, const double t, double dt,
		TracerStickerNames const& tracerstickernames, vector<Conserved3D> &extensives) const;

	double SuggestInverseTimeStep(void)const;

private:
	const Acceleration3D& acc_;
	const bool mass_flux_;
	mutable double dt_;
};

#endif // CONSFORCE3D_HPP