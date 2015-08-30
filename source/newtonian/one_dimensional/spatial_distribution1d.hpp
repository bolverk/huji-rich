/*! \file spatial_distribution1d.hpp
  \brief Abstract class for initial conditions
  \author Almog Yalinewich
 */

#ifndef SPATIAL_DISTRIBUTION_1D_HPP
#define SPATIAL_DISTRIBUTION_1D_HPP 1

//! \brief Base class for initial conditions
class SpatialDistribution1D
{
public:
  /*! \brief Calculates initial conditions
    \param x Position
    \result Value of the function at x
   */
  virtual double operator()(double x) const = 0;

  virtual ~SpatialDistribution1D(void);
};

//! \brief Uniform distribution
class Uniform: public SpatialDistribution1D
{
private:

  //! \brief Value at each point
  double _val;

public:
  /*! \brief Class constructor
    \param val Value at each point
   */
  explicit Uniform(double val);

  /*! \brief Calculates the value at each point
    \param x Position
    \return Value at x
   */
  double operator()(double x) const;
};

//! \brief Step distribution
class Step: public SpatialDistribution1D
{
private:

  //! \brief Value befor the step
  double _val1;

  //! \brief Value after the step
  double _val2;

  //! \brief Step position
  double _steppos;

public:
  /*! \brief Class constructor
    \param val1 Value befor the step
    \param val2 Value after the step
    \param steppos Position Step position
   */
  Step(double val1, double val2,
       double steppos);

  /*! \brief Evaluates the distribution at a certain spot
    \param x Position
    \return Value of the distribution at x
   */
  double operator()(double x) const;
};

//! \brief Two steps
class TwoSteps: public SpatialDistribution1D
{
private:
  //! \brief Value of the left step
  double vl;
  //! \brief Value of the right step
  double vr;
  //! \brief Value of the middle step
  double vm;
  //! \brief Positions of the left interface
  double lip;
  //! \brief Position of the right interface
  double rip;
public:
  /*! \brief Class constructor
    \param ivl Value on the left step
    \param ilip Position of the left interface
    \param ivm Value on the middle step
    \param irip Position of the right interface
    \param ivr Value on the right step
   */
  TwoSteps(double ivl, double ilip,
	   double ivm, double irip,
	   double ivr);
  /*! \brief Evaluates the distribution at a certain spot
    \param x Position
    \return Value of the distribution at x
   */
  double operator()(double x) const;
};

#endif // SPATIAL_DISTRIBUTION_1D_HPP
