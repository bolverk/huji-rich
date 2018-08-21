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

  /*! \brief Copy constructor
     \param origin Original
   */
  UniversalError(const UniversalError& origin);

  /*! \brief Appends std::string to the error message
    \param msg Message to append
   */
  void Append2ErrorMessage(std::string const& msg);

  /*! \brief Adds an entry to the list
    \param field Entry title
    \param value Value of entry
   */
  void AddEntry(std::string const& field,
		double value);

  /*! \brief Returns the error message
    \return Error message
   */
  std::string const& GetErrorMessage(void) const;

  /*! \brief Returns entry fields
    \return Textual fields
   */
  std::vector<std::string> const& GetFields(void) const;

  /*! \brief Returns entry values
    \return Numerical values
   */
  std::vector<double> const& GetValues(void) const;

  ~UniversalError(void);

private:

  std::string err_msg_;

  std::vector<std::string> fields_;

  std::vector<double> values_;

  UniversalError& operator=(const UniversalError&);
};

#endif // UNIVERSAL_ERROR_HPP
