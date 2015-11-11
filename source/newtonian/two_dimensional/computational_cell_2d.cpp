#include "computational_cell_2d.hpp"
#include "../../misc/utils.hpp"

typedef boost::container::flat_map<std::string, double> TracerMap;

ComputationalCell::ComputationalCell(void) :
density(0), pressure(0), velocity(Vector2D()), tracers(),
stickers(){}

ComputationalCell::ComputationalCell(ComputationalCell const& other) :
density(other.density), pressure(other.pressure), velocity(other.velocity), tracers(other.tracers),
stickers(other.stickers){}

ComputationalCell& ComputationalCell::operator=(ComputationalCell const& other)
{
	density = other.density;
	pressure = other.pressure;
	velocity = other.velocity;
	tracers = other.tracers;
	stickers = other.stickers;
	return *this;
}

ComputationalCell& ComputationalCell::operator+=(ComputationalCell const& other)
{
	this->density += other.density;
	this->pressure += other.pressure;
	this->velocity += other.velocity;
	assert(this->tracers.size() == other.tracers.size());
	for (size_t j = 0; j < this->tracers.size(); ++j)
	{
		assert((this->tracers.begin() + j)->first == (other.tracers.begin()+j)->first);
		(this->tracers.begin() + j)->second += (other.tracers.begin() + j)->second;
	}
	return *this;
}

ComputationalCell& ComputationalCell::operator-=(ComputationalCell const& other)
{
	this->density -= other.density;
	this->pressure -= other.pressure;
	this->velocity -= other.velocity;
	assert(this->tracers.size() == other.tracers.size());
	for (size_t j = 0; j < this->tracers.size(); ++j)
	{
		assert((this->tracers.begin() + j)->first == (other.tracers.begin() + j)->first);
		(this->tracers.begin() + j)->second -= (other.tracers.begin() + j)->second;
	}
	return *this;
}

ComputationalCell& ComputationalCell::operator*=(double s)
{
	this->density *= s;
	this->pressure *= s;
	this->velocity *= s;
	for (size_t j = 0; j < this->tracers.size(); ++j)
	{
		(this->tracers.begin() + j)->second *= s;
	}
	return *this;
}

void ComputationalCellAddMult(ComputationalCell &res, ComputationalCell const& other, double scalar)
{
	res.density += other.density*scalar;
	res.pressure += other.pressure*scalar;
	res.velocity += other.velocity*scalar;
	assert(res.tracers.size() == other.tracers.size());
	for (size_t j = 0; j < res.tracers.size(); ++j)
	{
		assert((res.tracers.begin() + j)->first == (other.tracers.begin() + j)->first);
		(res.tracers.begin() + j)->second += (other.tracers.begin() + j)->second*scalar;
	}
}

ComputationalCell operator+(ComputationalCell const& p1, ComputationalCell const& p2)
{
  ComputationalCell res(p1);
  res += p2;
  return res;
}

ComputationalCell operator-(ComputationalCell const& p1, ComputationalCell const& p2)
{
	ComputationalCell res(p1);
	res -= p2;
	return res;
}

ComputationalCell operator/(ComputationalCell const& p, double s)
{
	ComputationalCell res(p);
	res.density /= s;
	res.pressure /= s;
	for (boost::container::flat_map<std::string, double>::iterator it = res.tracers.begin(); it != res.tracers.end(); ++it)
		it->second /= s;
	res.velocity = res.velocity / s;
	return res;
}

ComputationalCell operator*(ComputationalCell const& p, double s)
{
	ComputationalCell res(p);
	res.density *= s;
	res.pressure *= s;
	for (boost::container::flat_map<std::string, double>::iterator it = res.tracers.begin(); it != res.tracers.end(); ++it)
		it->second *= s;
	res.velocity = res.velocity * s;
	return res;
}

ComputationalCell operator*(double s, ComputationalCell const& p)
{
	return p*s;
}
