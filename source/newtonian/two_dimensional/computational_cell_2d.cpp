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
      assert((this->tracers.begin() + static_cast<int>(j))->first == (other.tracers.begin()+static_cast<int>(j))->first);
      (this->tracers.begin() + static_cast<int>(j))->second += (other.tracers.begin() + static_cast<int>(j))->second;
    }
  return *this;
}

ComputationalCell& ComputationalCell::operator-=(ComputationalCell const& other)
{
  this->density -= other.density;
  this->pressure -= other.pressure;
  this->velocity -= other.velocity;
  assert(this->tracers.size() == other.tracers.size());
  /*for (size_t j = 0; j < this->tracers.size(); ++j)
    {
      assert((this->tracers.begin() + static_cast<int>(j))->first == (other.tracers.begin() + static_cast<int>(j))->first);
      (this->tracers.begin() + static_cast<int>(j))->second -= (other.tracers.begin() + static_cast<int>(j))->second;
    }*/
  boost::container::flat_map<string, double>::const_iterator it2 = other.tracers.begin();
  for (boost::container::flat_map<string, double>::iterator it = this->tracers.begin(); it != this->tracers.end(); ++it, ++it2)
	  it->second -= it2->second;
  return *this;
}

ComputationalCell& ComputationalCell::operator*=(double s)
{
  this->density *= s;
  this->pressure *= s;
  this->velocity *= s;
  for (size_t j = 0; j < this->tracers.size(); ++j)
    {
      (this->tracers.begin() + static_cast<int>(j))->second *= s;
    }
  return *this;
}

void ComputationalCellAddMult(ComputationalCell &res, ComputationalCell const& other, double scalar)
{
  res.density += other.density*scalar;
  res.pressure += other.pressure*scalar;
  res.velocity += other.velocity*scalar;
  assert(res.tracers.size() == other.tracers.size());
  /*const size_t nloop = res.tracers.size();
  for (size_t j = 0; j < nloop; ++j)
    {
      assert((res.tracers.begin() + static_cast<int>(j))->first == (other.tracers.begin() + static_cast<int>(j))->first);
      (res.tracers.begin() + static_cast<int>(j))->second += (other.tracers.begin() + static_cast<int>(j))->second*scalar;
    }*/
  boost::container::flat_map<string, double>::const_iterator it2 = other.tracers.begin();
  for(boost::container::flat_map<string, double>::iterator it = res.tracers.begin(); it != res.tracers.end(); ++it, ++it2)
		  it->second += it2->second*scalar;
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

void ReplaceComputationalCell(ComputationalCell & cell, ComputationalCell const& other)
{
  cell.density = other.density;
  cell.pressure = other.pressure;
  cell.velocity = other.velocity;
  assert(cell.tracers.size() == other.tracers.size());
  assert(cell.stickers.size() == other.stickers.size());
  /*const size_t ntracer = cell.tracers.size();
  const size_t nsticker = cell.stickers.size();
  for (size_t j = 0; j < ntracer; ++j)
    {
      assert((cell.tracers.begin() + static_cast<int>(j))->first == (other.tracers.begin() + static_cast<int>(j))->first);
      (cell.tracers.begin() + static_cast<int>(j))->second = (other.tracers.begin() + static_cast<int>(j))->second;
    }
  for (size_t j = 0; j < nsticker; ++j)
    {
      (cell.stickers.begin() + static_cast<int>(j))->second = (other.stickers.begin() + static_cast<int>(j))->second;
    }*/
  boost::container::flat_map<string, double>::const_iterator it2 = other.tracers.begin();
  for (boost::container::flat_map<string, double>::iterator it = cell.tracers.begin(); it != cell.tracers.end(); ++it, ++it2)
	  it->second = it2->second;
  boost::container::flat_map<string, bool>::const_iterator it3 = other.stickers.begin();
  for (boost::container::flat_map<string, bool>::iterator it = cell.stickers.begin(); it != cell.stickers.end(); ++it, ++it3)
	  it->second = it3->second;
}


Slope::Slope(void) :xderivative(ComputationalCell()), yderivative(ComputationalCell()) {}

Slope::Slope(ComputationalCell const & x, ComputationalCell const & y) : xderivative(x), yderivative(y)
{}


#ifdef RICH_MPI
size_t ComputationalCell::getChunkSize(void) const
{
  return 4 + tracers.size() + stickers.size();
}

vector<double> ComputationalCell::serialize(void) const
{
  vector<double> res (getChunkSize());
  res.at(0) = density;
  res.at(1) = pressure;
  res.at(2) = velocity.x;
  res.at(3) = velocity.y;
  size_t counter = 4;
  for(boost::container::flat_map<string,double>::const_iterator 
	it=tracers.begin();
      it!=tracers.end(); ++it){
    res.at(counter) = it->second;
    ++counter;
  }
  for(boost::container::flat_map<string,bool>::const_iterator
	it=stickers.begin();
      it!=stickers.end();
      ++it){
    res.at(counter) = it->second ? 0 : 1;
    ++counter;
  }
  assert(counter==res.size());
  return res;
}

void ComputationalCell::unserialize
(const vector<double>& data)
{
  assert(data.size()==getChunkSize());
  density = data.at(0);
  pressure = data.at(1);
  velocity.x = data.at(2);
  velocity.y = data.at(3);
  size_t counter = 4;
  for(boost::container::flat_map<string,double>::iterator 
	it=tracers.begin();
      it!=tracers.end(); ++it){
    it->second = data.at(counter);
    ++counter;
  }
  for(boost::container::flat_map<string,bool>::iterator
	it=stickers.begin();
      it!=stickers.end();
      ++it){
    it->second = data.at(counter)<0.5;
    ++counter;
  }
  assert(data.size()==counter);
}

size_t Slope::getChunkSize(void) const
{
	return xderivative.getChunkSize() * 2;
}

vector<double> Slope::serialize(void) const
{
	vector<double> res(getChunkSize());
	vector<double> temp(xderivative.serialize());
	std::copy(temp.begin(), temp.end(), res.begin());
	temp = yderivative.serialize();
	std::copy(temp.begin(), temp.end(), res.begin() + xderivative.getChunkSize());
	return res;
}

void Slope::unserialize(const vector<double>& data)
{
	size_t size = xderivative.getChunkSize();
	xderivative.unserialize(vector<double>(data.begin(), data.begin() + size));
	yderivative.unserialize(vector<double>(data.begin() + size, data.end()));
}

#endif // RICH_MPI
