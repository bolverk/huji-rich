/*! \file cell_updater_3d.hpp
  \brief Abstract class for cell update scheme
  \author Almog Yalinewich
 */

#ifndef CELL_UPDATER_HPP
#define CELL_UPDATER_HPP 1

#include "computational_cell.hpp"
#include "conserved_3d.hpp"
#include "../../3D/GeometryCommon/Tessellation3D.hpp"
#include "../two_dimensional/computational_cell_2d.hpp"
#include "../common/equation_of_state.hpp"

//! \brief Abstract clas for cell update scheme
class CellUpdater3D
{
public:

	/*! \brief Calculates the computational cell
	\param res The new cells given as output
	  \param extensives The extensive conserved variables
	  \param eos Equation of state
	  \param tess The tessellation
	  \param tracerstickernames The names of the tracers and stickers
	  \return Computational cell
	 */
	virtual void operator() (vector<ComputationalCell3D> &res, EquationOfState const& eos,
		const Tessellation3D& tess,vector<Conserved3D>& extensives,TracerStickerNames const& tracerstickernames)const = 0;

		//! \brief Class destructor
		virtual ~CellUpdater3D(void);
};

/*!
\brief Calculates velocity from extensive in SR
\param cell The extensive variable
\param G The adiabatic index
\return The velocity
*/
double GetVelocity(Conserved3D const& cell, double G);
#endif // CELL_UPDATER_HPP
