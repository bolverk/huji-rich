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

class Rubric: public Index2FileName
{
public:

  Rubric(const string& prefix,
	 const string& postfix);

  string operator()(int index);

private:
  const string prefix_;
  const string postfix_;
};

#endif // INDEX2FILENAME_HPP
