#include "SineSepar.hpp"

SineSepar::SineSepar(double Amp,double f,double off,double above,double under):
Amplitude(Amp),freq(f),offset(off),above_(above),under_(under)
{}

SineSepar::~SineSepar(void)
{}

double SineSepar::operator()(Vector2D const& r) const
{
	if(r.y>offset+Amplitude*sin(freq*r.x))
		return above_;
	else
		return under_;
}
