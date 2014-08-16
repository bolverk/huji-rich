#include "multiple_diagnostics.hpp"

MultipleDiagnostics::MultipleDiagnostics(void):
  diag_list() {}

void MultipleDiagnostics::operator()(const hdsim& sim)
{
  for(size_t i=0;i<diag_list.size();++i)
    (*diag_list[i])(sim);
}
