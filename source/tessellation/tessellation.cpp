#include "tessellation.hpp"
#include "../misc/universal_error.hpp"
#include <cmath>

int Tessellation::GetOriginalIndex(int /*point*/) const
{
	throw UniversalError("Method is only availible for Periodic boundaries");
}

Tessellation::~Tessellation()
{}
