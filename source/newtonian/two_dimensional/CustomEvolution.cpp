#include "CustomEvolution.hpp"
#include "../../misc/universal_error.hpp"

bool CustomEvolution::ShouldForceTracerReset(void)const
{
	return true;
}

bool CustomEvolution::isRelevantToInterpolation(void) const
{
  return false;
}

Vector2D CustomEvolution::CalcVelocity(int /*index*/,
				Tessellation const& /*tessellation*/,
				vector<Primitive> const& /*primitives*/,double /*time*/)
{
	return Vector2D(0,0);
}

CustomEvolution::~CustomEvolution(void) {}

bool CustomEvolution::TimeStepRelevant(void)const
{
  return true;
}

CustomEvolutionManager::CustomEvolutionManager(void):
  vtable_(1,NULL) {}

namespace
{
  bool check_for_duplicates(const vector<CustomEvolution*>& vtable,
			    CustomEvolution* cep)
  {
    for(size_t i=0;i<vtable.size();++i){
      if(vtable[i]==cep)
	return true;
    }
    return false;
  }
}

void CustomEvolutionManager::addCustomEvolution(CustomEvolution* cep)
{
  assert(!check_for_duplicates(vtable_,cep) &&
	 "Error in CustomEvolutionManager::addCustomEvolution: method already exists");

  vtable_.push_back(cep);
}

size_t CustomEvolutionManager::getIndex(CustomEvolution* cep) const
{
  for(size_t i=0;i<vtable_.size();++i){
    if(vtable_[i]==cep)
      return i;
  }
  throw UniversalError("Error in CustomEvolutionManager::getIndex: Method absent from vtable");
}

CustomEvolution* CustomEvolutionManager::getFunction(size_t i) const
{
  return vtable_[i];
}
