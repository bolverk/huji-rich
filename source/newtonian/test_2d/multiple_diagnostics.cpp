#include "multiple_diagnostics.hpp"

MultipleDiagnostics::MultipleDiagnostics(const vector<DiagnosticFunction*>& diag_list):
  diag_list_(diag_list) {}

void MultipleDiagnostics::operator()(const hdsim& sim)
{
  for(size_t i=0;i<diag_list_.size();++i)
    (*diag_list_[i])(sim);
}

MultipleDiagnostics::~MultipleDiagnostics(void)
{
  for(size_t i=0;i<diag_list_.size();++i)
    delete diag_list_[i];
}
