#include "PowerLawDistribution.hpp"

PowerLawDistribution::PowerLawDistribution(double A,double b,Vector2D const& center)
	:A_(A),b_(b),center_(center){}

double PowerLawDistribution::operator()(Vector2D const& point) const
{
	return A_*pow(abs(point-center_),b_);
}
