#include "CustomEvolution.hpp"

CustomEvolution::~CustomEvolution(void) {}

bool CustomEvolution::flux_indifferent(void) const
{
  return false;
}

bool CustomEvolution::TimeStepRelevant(void)const
{
	return true;
}