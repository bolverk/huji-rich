#include "conserved_3d.hpp"

using std::size_t;

Conserved3D::Conserved3D(void) :
	mass(0), momentum(), energy(0), internal_energy(0), tracers() {}

Conserved3D::Conserved3D(double mass_i,
	const Vector3D& momentum_i,
	double energy_i, double internal_energy_i) :
	mass(mass_i), momentum(momentum_i), energy(energy_i), internal_energy(internal_energy_i), tracers() {}

Conserved3D::Conserved3D(double mass_i,
	const Vector3D& momentum_i,
	double energy_i, double internal_energy_i,
	const std::array<double, MAX_TRACERS >& tracers_i) :
	mass(mass_i), momentum(momentum_i),
	energy(energy_i), internal_energy(internal_energy_i), tracers(tracers_i) {}

namespace
{
	std::array<double, MAX_TRACERS> operator*(double s, const std::array<double, MAX_TRACERS>& v)
	{
		std::array<double, MAX_TRACERS> res;
		for (size_t i = 0; i < v.size(); ++i)
			res[i] = s * v[i];
		return res;
	}

  /*
	std::array<double, MAX_TRACERS> operator/(const std::array<double, MAX_TRACERS>& v, double s)
	{
		std::array<double, MAX_TRACERS> res;
		double s_1 = 1.0 / s;
		for (size_t i = 0; i < v.size(); ++i)
			res[i] = v[i] * s_1;
		return res;
	}
  */
}

Conserved3D& Conserved3D::operator-=(const Conserved3D& diff)
{
	mass -= diff.mass;
	momentum -= diff.momentum;
	energy -= diff.energy;
	internal_energy -= diff.internal_energy;
	for (size_t i = 0; i < MAX_TRACERS; ++i)
		tracers[i] -= diff.tracers[i];
	return *this;
}

Conserved3D& Conserved3D::operator+=(const Conserved3D& diff)
{
	mass += diff.mass;
	momentum += diff.momentum;
	energy += diff.energy;
	internal_energy += diff.internal_energy;
	for (size_t i = 0; i < tracers.size(); ++i)
		tracers[i] += diff.tracers[i];
	return *this;
}

#ifdef RICH_MPI
size_t Conserved3D::getChunkSize(void) const
{
	return 6 + tracers.size();
}

vector<double> Conserved3D::serialize(void) const
{
	vector<double> res(getChunkSize());
	res.at(0) = mass;
	res.at(1) = energy;
	res.at(2) = momentum.x;
	res.at(3) = momentum.y;
	res.at(4) = momentum.z;
	res.at(5) = internal_energy;
	size_t counter = 6;
	//size_t N = tracers.size();
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		res[j + counter] = tracers[j];
	return res;
}

void Conserved3D::unserialize(const vector<double>& data)
{
	assert(data.size() == getChunkSize());
	mass = data.at(0);
	energy = data.at(1);
	momentum.x = data.at(2);
	momentum.y = data.at(3);
	momentum.z = data.at(4);
	internal_energy = data.at(5);
	size_t counter = 6;
	//size_t N = tracers.size();
	for (size_t j = 0; j < MAX_TRACERS; ++j)
		tracers[j] = data.at(counter + j);
}
#endif

Conserved3D operator*(double s, const Conserved3D& c)
{
	return Conserved3D(s*c.mass,
		s*c.momentum,
		s*c.energy, s*c.internal_energy,
		s*c.tracers);
}

Conserved3D operator*(const Conserved3D& c, double s)
{
	return Conserved3D(s*c.mass,
		s*c.momentum,
		s*c.energy, s*c.internal_energy,
		s*c.tracers);
}


Conserved3D operator/(const Conserved3D& c, double s)
{
	double s_1 = 1.0 / s;
	return Conserved3D(c.mass * s_1,
		c.momentum * s_1,
		c.energy * s_1, c.internal_energy * s_1,
		s_1 * c.tracers);
}

void PrimitiveToConserved(ComputationalCell3D const& cell, double vol, Conserved3D &res)
{
	res.mass = cell.density*vol;
	res.momentum = cell.velocity;
	res.momentum *= res.mass;
	res.internal_energy = res.mass*cell.internal_energy;
	res.energy = res.mass*0.5*ScalarProd(cell.velocity, cell.velocity) + res.internal_energy;
	//size_t N = cell.tracers.size();
	//res.tracers.resize(N);
	for (size_t i = 0; i < MAX_TRACERS; ++i)
		res.tracers[i] = cell.tracers[i] * res.mass;
}

void PrimitiveToConservedSR(ComputationalCell3D const& cell, double vol, Conserved3D &res, EquationOfState const& eos, TracerStickerNames const& tsn)
{
	double gamma = 1 / std::sqrt(1 - ScalarProd(cell.velocity, cell.velocity));
	res.mass = cell.density*vol*gamma;
	const double enthalpy = eos.dp2e(cell.density, cell.pressure, cell.tracers, tsn.tracer_names);
	res.internal_energy = enthalpy * res.mass;
	if (fastabs(cell.velocity) < 1e-5)
		res.energy = (gamma*enthalpy + 0.5*ScalarProd(cell.velocity, cell.velocity))* res.mass - cell.pressure*vol;
	else
		res.energy = (gamma*enthalpy + (gamma - 1))* res.mass - cell.pressure*vol;
	res.momentum = res.mass * (enthalpy + 1)*gamma*cell.velocity;
	size_t N = cell.tracers.size();
	//res.tracers.resize(N);
	for (size_t i = 0; i < N; ++i)
		res.tracers[i] = cell.tracers[i] * res.mass;
}

Conserved3D operator+(Conserved3D const& p1, Conserved3D const& p2)
{
	Conserved3D res(p1);
	res += p2;
	return res;
}

Conserved3D operator-(Conserved3D const& p1, Conserved3D const& p2)
{
	Conserved3D res(p1);
	res -= p2;
	return res;
}

Conserved3D& Conserved3D::operator*=(double s)
{
	this->mass *= s;
	this->momentum *= s;
	this->energy *= s;
	this->internal_energy *= s;
	size_t N = this->tracers.size();
	for (size_t j = 0; j < N; ++j)
		this->tracers[j] *= s;
	return *this;
}
