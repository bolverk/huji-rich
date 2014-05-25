#include <stdlib.h>
#include <cmath>
#include <vector>
#include <boost/array.hpp>
using namespace std;
using namespace boost;

/*
 * Returns the absolute value of the number given.
 */
double absolute(double const& a);

/*
 * Calculates the sum of a and b.
 * res - the computed result, err - the roundoff error.
 */
void fastTwoSum(double a, double b, double& res, double& err);

/*
 * Subtracts b from a.
 * res - the computed result, err - the roundoff error.
 */
void fastTwoDiff(double a, double b, double& res, double& err);

/*
 * Calculates the sum of a and b.
 * res - the computed result, err - the roundoff error.
 */
void twoSum(double a, double b, double& res, double& err);

/*
 * Subtracts b from a.
 * res - the computed result, err - the roundoff error.
 */
void twoDiff(double a, double b, double& res, double& err);

/*
 * Splits a given number into two, Used for multiplication.
 */
void split(double num, double& high, double& low);

/*
 * Multiplies a and b.
 * res - the computed result, err - the roundoff error.
 */
void twoProduct(double a, double b, double& res, double& err);

/*
 * Calculates the square of a number.
 * res - the computed result, err - the roundoff error.
 */
void square(double num, double& res, double& err);

/*
 * Calculates the sum of a two-expansion and a double.
 */
boost::array<double,3> twoOneSum(boost::array<double,2> const& a, double b);

/*
 * Calculates the difference between a two-expansion and a double.
 */
boost::array<double,3> twoOneDiff(boost::array<double,2> const& a, double b);

/*
 * Calculates the sum of two two-expansions.
 */
vector<double> twoTwoSum(boost::array<double,2> const& a,boost::array<double,2> const& b);

/*
 * Calculates the difference between two two-expansions.
 */
vector<double> twoTwoDiff(boost::array<double,2> const& a, boost::array<double,2> const& b);

/*
 * Adds a scalar to an existing expansion.
 */
vector<double> growExpansionZeroElim(vector<double> const& e, double b);

/*
 * Adds up two expansions.
 */
vector<double> expansionSumZeroElim(vector<double> const& e, vector<double> const& f);

/*
 * Adds up two expansions.
 */
vector<double> fastExpansionSumZeroElim(vector<double> const& e, vector<double> const& f);

/*
 * Adds up two expansions.
 */
vector<double> linearExpansionSumZeroElim(vector<double> const& e, vector<double> const& f);

/*
 * Multiplies a scalar by an expansion.
 */
void scaleExpansionZeroElim(vector<double> const& e, double b,
	vector<double> &result);

/*
 * Compresses an expansion.
 */
vector<double> compress(vector<double> const& e);

/*
 * Calculate a double precision approximation of the expansion.
 */
double estimate(vector<double> const& e);
