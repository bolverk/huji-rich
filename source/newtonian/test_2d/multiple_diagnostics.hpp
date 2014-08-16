/*! \file multiple_diagnostics.hpp
  \author Almog Yalinewich
  \brief Enables using more than a single diagnostic function
 */

#ifndef MULTIPLE_DIAGNOSTICS_HPP
#define MULTIPLE_DIAGNOSTICS_HPP 1

#include "main_loop_2d.hpp"

class MultipleDiagnostics: public DiagnosticFunction
{
public:

  vector<DiagnosticFunction*> diag_list;

  MultipleDiagnostics(void);

  void operator()(const hdsim& sim);
};

#endif // MULTIPLE_DIAGNOSTICS_HPP
