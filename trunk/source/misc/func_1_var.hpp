/*! \file func_1_var.hpp
  \brief Class that represents a mathematical function of a single variable
  \author Almog Yalinewich
*/
#ifndef FUNC_1_VAR_HPP
#define FUNC_1_VAR_HPP 1

//! \brief Scalar function of a single variable
class Func1Var
{
public:

  /*! \brief Evaluates the function
    \param x Independent variable
	\return The evaluated value
   */
  virtual double eval(double x) const = 0;
  //! \brief virtual destructor
  virtual ~Func1Var(void);
};

#endif // FUNC_1_VAR_HPP
