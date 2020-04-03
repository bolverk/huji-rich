#include <cmath>
#include "cell_updater_1d.hpp"

CellUpdater1D::~CellUpdater1D(void) {}

SimpleCellUpdater1D::SimpleCellUpdater1D(void) {}

namespace{
    Primitive retrive_single_cell
    (const Conserved& intensive,
     const EquationOfState& eos)
     {
         const double density = intensive.Mass;
         const Vector2D velocity = intensive.Momentum/intensive.Mass;
         const double kinetic_energy = 0.5*pow(abs(velocity),2);
         const double total_energy = intensive.Energy/intensive.Mass;
         const double thermal_energy = total_energy - kinetic_energy;
         const double pressure = eos.de2p(density, thermal_energy);
         const double sound_speed = eos.dp2c(density, pressure);
         return Primitive(density, 
                          pressure,
                          velocity,
                          thermal_energy,
                          sound_speed);
     }
}

vector<Primitive> SimpleCellUpdater1D::operator()
(const vector<Conserved>& intensives,
const vector<Conserved>& /*extensives*/,
const vector<ComputationalCell>& /*old*/,
const EquationOfState& eos) const
{
    vector<Primitive> res(intensives.size());
    for(size_t i=0;i<res.size();++i)
        res.at(i) = retrive_single_cell(intensives.at(i), eos);
    return res;
}

SimpleCellUpdater1D::~SimpleCellUpdater1D(void) {}
