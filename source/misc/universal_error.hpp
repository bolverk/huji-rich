/*! \file universal_error.hpp
  \brief A class for storing error and debug information
  \author Almog Yalinewich
 */

#ifndef UNIVERSAL_ERROR_HPP
#define UNIVERSAL_ERROR_HPP 1

#include <string>
#include <vector>

/*! \brief Container for error reports
 */
class UniversalError
{
public:

  /*! \brief Class constructor
    \param err_msg Error message
   */
  explicit UniversalError(std::string const& err_msg);

  /*! \brief Appends std::string to the error message
    \param msg Message to append
   */
  void Append2ErrorMessage(std::string const& msg);

  /*! \brief Adds an entry to the list
   */
  void AddEntry(std::string const& field,
		double value);

  /*! \brief Returns the error message
   */
  std::string const& GetErrorMessage(void) const;

  /*! \brief Returns entry fields
   */
  std::vector<std::string> const& GetFields(void) const;

  /*! \brief Returns entry values
   */
  std::vector<double> const& GetValues(void) const;

  ~UniversalError(void);

private:

  std::string err_msg_;

  std::vector<std::string> fields_;

  std::vector<double> values_;
};

/*! \brief Prints the contents of the error
\param eo The error object
*/
void DisplayError(UniversalError const& eo);

#endif // UNIVERSAL_ERROR_HPP
