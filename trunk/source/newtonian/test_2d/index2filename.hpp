/*! \brief index2filename.hpp
  \brief Converts an index to a filename
  \author Almog Yalinewich
 */

#ifndef INDEX2FILENAME_HPP
#define INDEX2FILENAME_HPP 1

#include <string>

using std::string;

//! \brief Class for generating sequential file names
class Index2FileName
{
public:

  /*! \brief Generate a file name
    \param index Index
    \return File name
   */
  virtual string operator()(int index) = 0;

  //! \brief class destructor
  virtual ~Index2FileName(void);
};

//! \brief Class for creating files names using the pattern prefix+index+postfix
class Rubric: public Index2FileName
{
public:

  /*! \brief Class constructor
    \param prefix Prefix
    \param postfix Postfix
   */
  Rubric(const string& prefix,
	 const string& postfix);

  /*! \brief Constructs a file name
    \param index File index
    \return File name
   */
  string operator()(int index);

private:
  const string prefix_;
  const string postfix_;
};

#endif // INDEX2FILENAME_HPP
