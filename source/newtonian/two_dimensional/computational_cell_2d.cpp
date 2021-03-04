#include "computational_cell_2d.hpp"
#include "../../misc/utils.hpp"

using namespace std;

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
  size_t N = this->tracers.size();
  for (size_t j = 0; j < N; ++j)
    this->tracers[j] += other.tracers[j];
  return *this;
}

ComputationalCell& ComputationalCell::operator-=(ComputationalCell const& other)
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

ComputationalCell& ComputationalCell::operator*=(double s)
{
  this->density *= s;
  this->pressure *= s;
  this->velocity *= s;
  size_t N = this->tracers.size();
  for (size_t j = 0; j < N; ++j)
    this->tracers[j] *= s;
  return *this;
}

void ComputationalCellAddMult(ComputationalCell &res, ComputationalCell const& other, double scalar)
{
  res.density += other.density*scalar;
  res.pressure += other.pressure*scalar;
  res.velocity += other.velocity*scalar;
  assert(res.tracers.size() == other.tracers.size());
  size_t N = res.tracers.size();
  for (size_t j = 0; j < N; ++j)
    res.tracers[j] += other.tracers[j]*scalar;
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
  size_t N = res.tracers.size();
  for (size_t j = 0; j < N; ++j)
    res.tracers[j] /= s;
  res.velocity = res.velocity / s;
  return res;
}

ComputationalCell operator*(ComputationalCell const& p, double s)
{
  ComputationalCell res(p);
  res.density *= s;
  res.pressure *= s;
  size_t N = res.tracers.size();
  for (size_t j = 0; j < N; ++j)
    res.tracers[j] *= s;
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
  size_t N = cell.tracers.size();
  for (size_t j = 0; j < N; ++j)
    cell.tracers[j] = other.tracers[j];
  N = cell.stickers.size();
  for (size_t i = 0; i < N; ++i)
    cell.stickers[i] = other.stickers[i];
}


Slope::Slope(void) :xderivative(ComputationalCell()), yderivative(ComputationalCell()) {}

Slope::Slope(ComputationalCell const & x, ComputationalCell const & y) : xderivative(x), yderivative(y)
{}

TracerStickerNames::TracerStickerNames(void)
  : tracer_names(vector<string>()), sticker_names(vector<string>()) {}

TracerStickerNames& TracerStickerNames::operator=
(const TracerStickerNames& other)
{
  tracer_names = other.tracer_names;
  sticker_names = other.sticker_names;
  return *this;
}

TracerStickerNames::~TracerStickerNames(void) {}

TracerStickerNames::TracerStickerNames(TracerStickerNames const& other):tracer_names(other.tracer_names),sticker_names(other.sticker_names){}

TracerStickerNames::TracerStickerNames(std::vector<std::string> tracers, std::vector<std::string> stickers):
  tracer_names(tracers),sticker_names(stickers){}

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
  size_t N = tracers.size();
  for (size_t j = 0; j < N ; ++j)
    res[j+counter] = tracers[j];
  size_t N2 = stickers.size();
  for (size_t j = 0; j < N2; ++j)
    res[j + counter + N] = stickers[j] ? 1 : 0;
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
  size_t N = tracers.size();
  for (size_t j = 0; j < N; ++j)
    tracers[j] = data.at(counter + j);
  size_t N2 = stickers.size();
  for (size_t i = 0; i < N2; ++i)
    stickers[i] = data.at(counter + N + i)>0.5;
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
