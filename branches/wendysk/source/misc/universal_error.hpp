#ifndef UNIVERSAL_ERROR_HPP
#define UNIVERSAL_ERROR_HPP 1

#include <string>
#include <vector>

using namespace std;

/*! \brief Container for error reports
 */
class UniversalError
{
public:

  /*! \brief Class constructor
    \param err_msg Error message
   */
  UniversalError(string const& err_msg);

  /*! \brief Appends string to the error message
    \param msg Message to append
   */
  void Append2ErrorMessage(string const& msg);

  /*! \brief Adds an entry to the list
   */
  void AddEntry(string const& field,
		double value);

  /*! \brief Returns the error message
   */
  string const& GetErrorMessage(void) const;

  /*! \brief Returns entry fields
   */
  vector<string> const& GetFields(void) const;

  /*! \brief Returns entry values
   */
  vector<double> const& GetValues(void) const;

  ~UniversalError(void);
  
private:

  string err_msg_;

  vector<string> fields_;

  vector<double> values_;

};

#endif // UNIVERSAL_ERROR_HPP
