/*! \file multiple_diagnostics.hpp
  \author Almog Yalinewich
  \brief Enables using more than a single diagnostic function
 */

#ifndef MULTIPLE_DIAGNOSTICS_HPP
#define MULTIPLE_DIAGNOSTICS_HPP 1

#include "main_loop_2d.hpp"

//! \brief A class that enables using multiple diagnostics simultaneously
class MultipleDiagnostics: public DiagnosticFunction
{
public:

  //! \brief List of diagnostics
  const vector<DiagnosticFunction*> diag_list_;

  MultipleDiagnostics(const vector<DiagnosticFunction*>& diag_list);

  void operator()(const hdsim& sim);

  ~MultipleDiagnostics(void);
};

#endif // MULTIPLE_DIAGNOSTICS_HPP
