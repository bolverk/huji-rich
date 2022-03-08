/*! \file default_cell_updater.hpp
  \brief Default cell update scheme
  \author Almog Yalinewich
 */

#ifndef DEFAULT_CELL_UPDATER_HPP
#define DEFAULT_CELL_UPDATER_HPP 1

#include "cell_updater_3d.hpp"

 //! \brief Default scheme for cell update
class DefaultCellUpdater : public CellUpdater3D
{
public:

	//! \brief Class constructor
  //! \param SR Special relativity flag
  //! \param G Correction to adiabatic index
  //! \param includes_temperature Flag if to compute the temperature as well or not
	DefaultCellUpdater(bool SR = false,double G=0, bool const includes_temperature = false);

	void operator()(vector<ComputationalCell3D> &res, EquationOfState const& eos,
		const Tessellation3D& tess, vector<Conserved3D>& extensives) const override;
private:
	const bool SR_;
	const double G_;
	const bool includes_temperature_;
	mutable size_t entropy_index_;
};

#endif // DEFAULT_CELL_UPDATER_HPP
