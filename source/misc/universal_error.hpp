/*! \file universal_error.hpp
  \brief A class for storing error and debug information
  \author Almog Yalinewich
 */

#ifndef UNIVERSAL_ERROR_HPP
#define UNIVERSAL_ERROR_HPP 1

#include <string>
#include <vector>

using std::string;
using std::vector;
using std::pair;

/*! \brief Container for error reports
 */
class UniversalError
{
public:

  /*! \brief Class constructor
    \param err_msg Error message
   */
  explicit UniversalError(const string& err_msg);

  /*! \brief Appends std::string to the error message
    \param msg Message to append
   */
  void Append2ErrorMessage(std::string const& msg);

  /*! \brief Adds an entry to the list
    \param field Field name
    \param value value
   */
  void addEntry(std::string const& field,
		double value);

  /*! \brief Returns the error message
    \return Error message
   */
  std::string const& getErrorMessage(void) const;

  /*! \brief Returns entry fields and values
    \return List of fields and values
   */
  std::vector<pair<string, double> > const& getFields(void) const;

  ~UniversalError(void);

  /*! \brief Copy constructor
    \param eo Source
   */
  UniversalError(const UniversalError& eo);

private:

  string err_msg_;

  vector<pair<string, double> > fields_;
};

/*! \brief Prints the contents of the error
\param eo The error object
*/
void reportError(const UniversalError& eo);

#endif // UNIVERSAL_ERROR_HPP
