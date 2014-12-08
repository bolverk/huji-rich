/*! \file exactmath.hpp
  \brief Exact math
  \author Elad Steinberg
 */

#ifndef EXACTMATH_HPP
#define EXACTMATH_HPP 1

#include <stdlib.h>
#include <cmath>
#include <vector>
#include <boost/array.hpp>

using std::vector;
using std::min;

/*! \brief Calculates the sum of two numbers
  \param a First number
  \param b Second number
  \param res Output
  \param err Roundoff error
 */
void fastTwoSum(double a, double b, double& res, double& err);

/*! \brief Subtracts two numbers
   \param a First number
   \param b Second number
   \param res Output
   \param err Roundoff error
 */
void fastTwoDiff(double a, double b, double& res, double& err);

/*! \brief Calculates the sum of a and b.
   \param a First number
   \param b Second number
   \param res Result
   \param err Roundoff error
 */
void twoSum(double a, double b, double& res, double& err);

/*! \brief Difference between two numbers
   \param a First number
   \param b Second number
   \param res Result
   \param err Roundoff error
 */
void twoDiff(double a, double b, double& res, double& err);

/*! \brief Splits a given number into two, Used for multiplication.
   \param num Number
   \param high Higher part
   \param low Lower part
 */
void split(double num, double& high, double& low);

/*! \brief Product of two number
   \param a First number
   \param b Second number
   \param res Result
   \param err Error
 */
void twoProduct(double a, double b, double& res, double& err);

/*! \brief Calculates the square of a number.
   \param num Number
   \param res Result
   \param err Roundoff error
 */
void square(double num, double& res, double& err);

/*! \brief Calculates the sum of a two-expansion and a double.
   \param a Two expansion
   \param b A number
   \return Sum
 */
boost::array<double,3> twoOneSum(boost::array<double,2> const& a, double b);

/*! \brief Calculates the difference between a two-expansion and a double.
   \param a Two expansion
   \param b Number
   \return Difference
 */
boost::array<double,3> twoOneDiff(boost::array<double,2> const& a, double b);

/*! \brief Calculates the sum of two two-expansions.
   \param a First two expansion
   \param b Second two expansio
   \return sum
 */
vector<double> twoTwoSum(boost::array<double,2> const& a,boost::array<double,2> const& b);

/*! \brief Calculates the difference between two two-expansions.
   \param a First two expansion
   \param b Second two expansion
   \return Difference
 */
vector<double> twoTwoDiff(boost::array<double,2> const& a, boost::array<double,2> const& b);

/*! \brief Adds a scalar to an existing expansion.
  \param e Expansion
  \param b Scalar
  \return Expansion
 */
vector<double> growExpansionZeroElim(vector<double> const& e, double b);

/*! \brief Adds up two expansions.
   \param e First expansion
   \param f Second expansion
   \return Sum
 */
vector<double> expansionSumZeroElim(vector<double> const& e, vector<double> const& f);

/*! \brief Adds up two expansions.
  \param e First expansion
  \param f Second expansion
  \return Sum
 */
vector<double> fastExpansionSumZeroElim(vector<double> const& e, vector<double> const& f);

/*! \brief Adds up two expansions.
  \param e First expansion
  \param f Second expansion
  \return Expansion
 */
vector<double> linearExpansionSumZeroElim(vector<double> const& e, vector<double> const& f);

/*! \brief Multiplies a scalar by an expansion.
  \param e Expansion
  \param b Scalar
  \param result Result
 */
void scaleExpansionZeroElim(vector<double> const& e, double b,
	vector<double> &result);

/*! \brief Compresses an expansion.
   \param e Expansion
   \return Expansion
 */
vector<double> compress(vector<double> const& e);

/*! \brief Calculate a double precision approximation of the expansion.
  \param e Expansion
  \return A number
 */
double estimate(vector<double> const& e);

#endif //EXACTMATH_HPP
