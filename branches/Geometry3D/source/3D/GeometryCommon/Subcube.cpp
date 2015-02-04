//\file Subcube.cpp
//\brief Implementation of the subcube class - representing one of the 27 subcubes
//\author Itay Zandbank

#include "Subcube.hpp"
#include <set>
using namespace std;

Subcube::Subcube(const char offsets[3])
{
	for (int i = 0; i < 3; i++)
	{
		if (offsets[i] != MINUS && offsets[i] != CENTER && offsets[i] != PLUS)
			throw invalid_argument("Offset must be either '-', ' ' or '+'");

		_offsets[i] = offsets[i];
	}
}

bool operator<(const Subcube &sc1, const Subcube &sc2)
{
	return sc1.Num() < sc2.Num();
}

bool operator==(const Subcube &sc1, const Subcube &sc2)
{
	return sc1.Num() == sc2.Num();
}


set<Subcube> Subcube::_all;
const set<Subcube> &Subcube::all()
{
	if (_all.empty())
	{
		const char *offsets[] = {
			"---", "-- ", "--+", "- -", "-  ", "- +", "-+-", "-+ ", "-++",
			" --", " - ", " -+", "  -", "  +", " +-", " + ", " ++",
			"+--", "+- ", "+-+", "+ -", "+  ", "+ +", "++-", "++ ", "+++",
		};
		for (int i = 0; i < 26; i++)
			_all.insert(Subcube(offsets[i]));
	}

	return _all;
}