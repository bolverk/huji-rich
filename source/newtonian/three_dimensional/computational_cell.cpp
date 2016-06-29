#include "computational_cell.hpp"

ComputationalCell3D::ComputationalCell3D(void):
  density(0), pressure(0), velocity(), tracers(),stickers() {}

ComputationalCell3D::ComputationalCell3D(double density_i,
				     double pressure_i,
				     const Vector3D& velocity_i):
  density(density_i), pressure(pressure_i), 
  velocity(velocity_i), tracers(),stickers() {}

ComputationalCell3D::ComputationalCell3D(double density_i,
				     double pressure_i,
				     const Vector3D& velocity_i,
				     const vector<double>& tracers_i,
					 const vector<bool>& stickers_i):
  density(density_i), pressure(pressure_i),
  velocity(velocity_i), tracers(tracers_i),stickers(stickers_i) {}


ComputationalCell3D& ComputationalCell3D::operator=(ComputationalCell3D const& other)
{
	density = other.density;
	pressure = other.pressure;
	velocity = other.velocity;
	tracers = other.tracers;
	stickers = other.stickers;
	return *this;
}

ComputationalCell3D& ComputationalCell3D::operator+=(ComputationalCell3D const& other)
{
	this->density += other.density;
	this->pressure += other.pressure;
	this->velocity += other.velocity;
	assert(this->tracers.size() == other.tracers.size());
	size_t N = this->tracers.size();
	for (size_t j = 0; j < N; ++j)
		this->tracers[j] += other.tracers[j];
	return *this;
}

ComputationalCell3D& ComputationalCell3D::operator-=(ComputationalCell3D const& other)
{
	this->density -= other.density;
	this->pressure -= other.pressure;
	this->velocity -= other.velocity;
	assert(this->tracers.size() == other.tracers.size());
	size_t N = this->tracers.size();
	for (size_t j = 0; j < N; ++j)
		this->tracers[j] -= other.tracers[j];
	return *this;
}

ComputationalCell3D& ComputationalCell3D::operator*=(double s)
{
	this->density *= s;
	this->pressure *= s;
	this->velocity *= s;
	size_t N = this->tracers.size();
	for (size_t j = 0; j < N; ++j)
		this->tracers[j] *= s;
	return *this;
}

void ComputationalCellAddMult(ComputationalCell3D &res, ComputationalCell3D const& other, double scalar)
{
	res.density += other.density*scalar;
	res.pressure += other.pressure*scalar;
	res.velocity += other.velocity*scalar;
	assert(res.tracers.size() == other.tracers.size());
	size_t N = res.tracers.size();
	for (size_t j = 0; j < N; ++j)
		res.tracers[j] += other.tracers[j] * scalar;
}

ComputationalCell3D operator+(ComputationalCell3D const& p1, ComputationalCell3D const& p2)
{
	ComputationalCell3D res(p1);
	res += p2;
	return res;
}

ComputationalCell3D operator-(ComputationalCell3D const& p1, ComputationalCell3D const& p2)
{
	ComputationalCell3D res(p1);
	res -= p2;
	return res;
}

ComputationalCell3D operator/(ComputationalCell3D const& p, double s)
{
	ComputationalCell3D res(p);
	res.density /= s;
	res.pressure /= s;
	size_t N = res.tracers.size();
	for (size_t j = 0; j < N; ++j)
		res.tracers[j] /= s;
	res.velocity = res.velocity / s;
	return res;
}

ComputationalCell3D operator*(ComputationalCell3D const& p, double s)
{
	ComputationalCell3D res(p);
	res.density *= s;
	res.pressure *= s;
	size_t N = res.tracers.size();
	for (size_t j = 0; j < N; ++j)
		res.tracers[j] *= s;
	res.velocity = res.velocity * s;
	return res;
}

ComputationalCell3D operator*(double s, ComputationalCell3D const& p)
{
	return p*s;
}

void ReplaceComputationalCell(ComputationalCell3D & cell, ComputationalCell3D const& other)
{
	cell.density = other.density;
	cell.pressure = other.pressure;
	cell.velocity = other.velocity;
	assert(cell.tracers.size() == other.tracers.size());
	assert(cell.stickers.size() == other.stickers.size());
	size_t N = cell.tracers.size();
	for (size_t j = 0; j < N; ++j)
		cell.tracers[j] = other.tracers[j];
	N = cell.stickers.size();
	for (size_t i = 0; i < N; ++i)
		cell.stickers[i] = other.stickers[i];
}



Slope3D::Slope3D(void) : xderivative(ComputationalCell3D()), yderivative(ComputationalCell3D()), zderivative(ComputationalCell3D()) {}

Slope3D::Slope3D(ComputationalCell3D const & x, ComputationalCell3D const & y,ComputationalCell3D const & z) : xderivative(x), yderivative(y),zderivative(z)
{}
