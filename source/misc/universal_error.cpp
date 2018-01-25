#include "universal_error.hpp"
#include <iostream>

using namespace std;

UniversalError::UniversalError(string const& err_msg):
  err_msg_(err_msg),
  fields_(vector<string>()),
  values_(vector<double>()) {}

void UniversalError::Append2ErrorMessage(string const& msg)
{
  err_msg_ += msg;
}

void UniversalError::AddEntry(string const& field,
	      double value)
{
  fields_.push_back(field);
  values_.push_back(value);
}

string const& UniversalError::GetErrorMessage(void) const
{
  return err_msg_;
}

vector<string> const& UniversalError::GetFields(void) const
{
  return fields_;
}

vector<double> const& UniversalError::GetValues(void) const
{
  return values_;
}

UniversalError::~UniversalError(void) {}

void DisplayError(UniversalError const& eo)
{
	std::cout.precision(14);
	std::cout << eo.GetErrorMessage() << std::endl;
	for (size_t i = 0; i<eo.GetFields().size(); ++i)
		std::cout << eo.GetFields()[i] << " = " << eo.GetValues()[i] << std::endl;
}
