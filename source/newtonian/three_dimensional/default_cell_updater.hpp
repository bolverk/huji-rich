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
	DefaultCellUpdater(void);

	void operator()(vector<ComputationalCell3D> &res, EquationOfState const& eos,
		const Tessellation3D& tess, vector<Conserved3D>& extensives, 
		TracerStickerNames const& tracerstickernames) const;
private:
	mutable size_t entropy_index_;
};

#endif // DEFAULT_CELL_UPDATER_HPP
