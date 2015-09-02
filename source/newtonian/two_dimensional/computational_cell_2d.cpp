#include "computational_cell_2d.hpp"

ComputationalCell::ComputationalCell(void) :
density(0), pressure(0), velocity(Vector2D()), tracers(),
stickers(){}

ComputationalCell::ComputationalCell(ComputationalCell const& other) :
density(other.density), pressure(other.pressure), velocity(other.velocity), tracers(other.tracers),
stickers(other.stickers){}

ComputationalCell& ComputationalCell::operator+=(ComputationalCell const& other)
{
	this->density += other.density;
	this->pressure += other.pressure;
	this->velocity += other.velocity;
	for (boost::container::flat_map<std::string, double>::iterator it = this->tracers.begin();
		it != this->tracers.end(); ++it)
		it->second += other.tracers.find(it->first)->second;
	return *this;
}

ComputationalCell operator+(ComputationalCell const& p1, ComputationalCell const& p2)
{
	ComputationalCell res(p1);
	res.density += p2.density;
	res.pressure += p2.pressure;
	for (boost::container::flat_map<std::string, double>::iterator it = res.tracers.begin(); it != res.tracers.end(); ++it)
		it->second += p2.tracers.find(it->first)->second;
	res.velocity += p2.velocity;
	return res;
}

ComputationalCell operator-(ComputationalCell const& p1, ComputationalCell const& p2)
{
	ComputationalCell res(p1);
	res.density -= p2.density;
	res.pressure -= p2.pressure;
	for (boost::container::flat_map<std::string, double>::iterator it = res.tracers.begin(); it != res.tracers.end(); ++it)
		it->second -= p2.tracers.find(it->first)->second;
	res.velocity -= p2.velocity;
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
