#include "universal_error.hpp"
#include <iostream>
#include <algorithm>

using namespace std;

UniversalError::UniversalError(string const& err_msg):
  err_msg_(err_msg),
  fields_() {}

void UniversalError::Append2ErrorMessage(string const& msg)
{
  err_msg_ += msg;
}

void UniversalError::addEntry
(const string& field,
 double value)
{
  fields_.push_back({field, value});
}

const string& UniversalError::getErrorMessage(void) const
{
  return err_msg_;
}

const vector<pair<string, double> >& UniversalError::getFields(void) const
{
  return fields_;
}

UniversalError::~UniversalError(void) {}

UniversalError::UniversalError(const UniversalError& eo):
  err_msg_(eo.getErrorMessage()),
  fields_(eo.getFields()) {}

void reportError(UniversalError const& eo)
{
  std::cout.precision(14);
  std::cout << eo.getErrorMessage() << std::endl;
  for_each
    (eo.getFields().begin(),
     eo.getFields().end(),
     [](const pair<string, double>& f)
     {cout << f.first << " " << f.second << endl;});
}
