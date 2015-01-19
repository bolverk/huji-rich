/*! \brief index2filename.hpp
  \brief Converts an index to a filename
  \author Almog Yalinewich
 */

#ifndef INDEX2FILENAME_HPP
#define INDEX2FILENAME_HPP 1

#include <string>

using std::string;

class Index2FileName
{
public:

  virtual string operator()(int index) = 0;

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
