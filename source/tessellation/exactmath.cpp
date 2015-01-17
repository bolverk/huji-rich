#include "exactmath.hpp"
double const splitter = 134217729;

namespace {
  void fastTwoSumTail(double const a, double const b, double const x, double& y)
  {
    //    const double bVirt = x - a;
    y = b - (x - a);
  }
}

void fastTwoSum(double const a, double const b, double& res, double& err)
{
  res = a + b;
  fastTwoSumTail(a, b, res, err);
}

namespace {
  void fastTwoDiffTail(double const a, double const b, double const x, double& y)
  {
    const double bVirt = a - x;
    y = bVirt - b;
  }
}

void fastTwoDiff(double a, double b, double& res, double& err)
{
  res = a - b;
  fastTwoDiffTail(a, b, res, err);
}

namespace {
  double twoSumTail(double const a, double const b, double const x)
  {
    double bVirt = x - a;
    double aVirt = x - bVirt;
    bVirt = b - bVirt;
    aVirt = a - aVirt;
    return aVirt + bVirt;
  }
}

void twoSum(double a, double b, double& res, double& err)
{
  res = a + b;
  err = twoSumTail(a, b, res);
}

namespace {
  void twoDiffTail(double const a, double const b, double const x, double& y)
  {
    double bVirt = a - x;
    double aVirt = x + bVirt;
    bVirt -= b;
    aVirt = a - aVirt;
    y = aVirt + bVirt;
  }
}

void twoDiff(double a, double b, double& res, double& err)
{
  res = a - b;
  twoDiffTail(a, b, res, err);
}

void split(double num, double& high, double& low)
{
  double c = splitter * num;
  double temp = c - num;
  high = c - temp;
  low = num - high;
}

namespace {
  void twoProductTail(double const a, double const b, double const x, double& y)
  {
    double ahi, alo, bhi, blo;
    split(a, ahi, alo);
    split(b, bhi, blo);
    y = x - (ahi * bhi);
    y -= (alo * bhi);
    y -= (ahi * blo);
    y = (alo * blo) - y;
  }
}

void twoProduct(double a, double b, double& res, double& err)
{
  res = a * b;
  twoProductTail(a, b, res, err);
}

namespace {
  void squareTail(double const a, double const x, double& y)
  {
    double high, low;
    split(a, high, low);
    y = x - (high * high);
    y -= ((high + high) * low);
    y = (low * low) - y;
  }
}

void square(double num, double& res, double& err)
{
  res = num * num;
  squareTail(num, res, err);
}

boost::array<double,3> twoOneSum(boost::array<double,2> const& a, double b)
{
  boost::array<double,3> res;
  double tmp;
  twoSum(a[0], b, tmp, res[0]);
  twoSum(a[1], tmp, res[2], res[1]);
  return res;
}

boost::array<double,3> twoOneDiff(boost::array<double,2> const& a, double b)
{
  boost::array<double,3> res;
  double tmp;
  twoDiff(a[0], b, tmp, res[0]);
  twoSum(a[1], tmp, res[2], res[1]);
  return res;
}

vector<double> twoTwoSum(boost::array<double,2> const& a, boost::array<double,2> const& b)
{
  vector<double> res(4);
  double tmp1, tmp2, tmp3;
  twoSum(a[0], b[0], tmp1, res[0]);
  twoSum(a[1], tmp1, tmp2, tmp3);
  twoSum(tmp3, b[1], tmp1, res[1]);
  twoSum(tmp2, tmp1, res[3], res[2]);
  return res;
}

vector<double> twoTwoDiff(boost::array<double,2> const& a, boost::array<double,2> const& b)
{
  vector<double> res(4);
  double tmp1, tmp2, tmp3;
  twoDiff(a[1], b[1], tmp1, res[0]);
  twoSum(a[0], tmp1, tmp2, tmp3);
  twoDiff(tmp3, b[0], tmp1, res[1]);
  twoSum(tmp2, tmp1, res[3], res[2]);
  return res;
}

vector<double> growExpansionZeroElim(vector<double> const& e, double b)
{
  vector<double> ans;
  ans.reserve(8);
  double R1, R2, res; //Temp variables to store values.

  R1 = b; //The remainder to add to the expansion.
  for (size_t i = 0; i < e.size(); ++i)
    {
      twoSum(R1, e[i], R2, res); //Adding the remainder to the current number in the expansion.
      R1 = R2; //Updating remainder.
      if (res != 0.0)
	{ //Ignoring 0.
	  ans.push_back(res); //Creating the result expansion.
	}
    }
  if ((R1 != 0.0) || ans.empty())
    { //No numbers were previously added to result expansion.
      ans.push_back(R1);
    }
  return ans; //Returns the result vector.
}

vector<double> expansionSumZeroElim(vector<double> const& e, vector<double> const& f)
{
  vector<double> temp, ans;
  temp.reserve(8);
  ans.reserve(8);
  double Q, Qnew, current;

  Q = f[0];
  for (size_t i = 0; i < e.size(); i++)
    {
      twoSum(Q, e[i], Qnew, current);
      temp.push_back(current);
      Q = Qnew;
    }
  temp.push_back(Q);
  for (size_t i = 1; i < f.size(); i++)
    {
      Q = f[i];
      for (size_t j = i; j < temp.size(); j++)
	{
	  current = temp[j];
	  twoSum(Q, current, Qnew, temp[j]);
	  Q = Qnew;
	}
      temp.push_back(Q);
    }
  for (size_t i = 0; i < temp.size(); i++)
    {
      if (fabs(temp[i]) > 0)
	{
	  ans.push_back(temp[i]);
	}
    }
  if (ans.empty())
    {
      ans.push_back(0.0);
    }
  return ans;
}

vector<double> fastExpansionSumZeroElim(vector<double> const& e, vector<double> const& f)
{
  vector<double> ans;
  ans.reserve(10);
  double Q, Qnew;
  double hh;
  double eNow, fNow;
  unsigned int eIndex, fIndex;

  eNow = e[0];
  fNow = f[0];
  eIndex = fIndex = 0;
  if ((fNow > eNow) == (fNow > -eNow))
    {
      Q = eNow;
      eIndex++;
      if (eIndex <e.size())
	{
	  eNow = e[eIndex];
	}
    }
  else
    {
      Q = fNow;
      fIndex++;
      if (fIndex <f.size())
	{
	  fNow = f[fIndex];
	}
    }
  if ((eIndex < e.size()) && (fIndex < f.size()))
    {
      if ((fNow > eNow) == (fNow > -eNow))
	{
	  fastTwoSum(eNow, Q, Qnew, hh);
	  eIndex++;
	  if (eIndex <e.size())
	    {
	      eNow = e[eIndex];
	    }
	}
      else
	{
	  fastTwoSum(fNow, Q, Qnew, hh);
	  fIndex++;
	  if (fIndex <f.size())
	    {
	      fNow = f[fIndex];
	    }
	}
      Q = Qnew;
      if (hh != 0.0)
	{
	  ans.push_back(hh);
	}
      while ((eIndex < e.size()) && (fIndex <f.size()))
	{
	  if ((fNow > eNow) == (fNow > -eNow))
	    {
	      twoSum(Q, eNow, Qnew, hh);
	      eIndex++;
	      if (eIndex < e.size())
		{
		  eNow = e[eIndex];
		}
	    }
	  else
	    {
	      twoSum(Q, fNow, Qnew, hh);
	      fIndex++;
	      if (fIndex < f.size())
		{
		  fNow = f[fIndex];
		}
	    }
	  Q = Qnew;
	  if (hh != 0.0)
	    {
	      ans.push_back(hh);
	    }
	}
    }
  while (eIndex < e.size())
    {
      twoSum(Q, eNow, Qnew, hh);
      eIndex++;
      if (eIndex < e.size())
	{
	  eNow = e[eIndex];
	}
      Q = Qnew;
      if (hh != 0.0)
	{
	  ans.push_back(hh);
	}
    }
  while (fIndex < f.size())
    {
      twoSum(Q, fNow, Qnew, hh);
      fIndex++;
      if (fIndex < f.size())
	{
	  fNow = f[fIndex];
	}
      Q = Qnew;
      if (hh != 0.0)
	{
	  ans.push_back(hh);
	}
    }
  if ((Q != 0.0) || ans.empty())
    {
      ans.push_back(Q);
    }
  return ans;
}

vector<double> linearExpansionSumZeroElim(vector<double> const& e, vector<double> const& f)
{
  vector<double> ans;
  ans.reserve(8);
  double Q, q, hh;
  double Qnew;
  double R;
  double g0;

  double enow = e.front();
  double fnow = f.front();
  size_t eindex = 0;
  size_t findex = 0;
  if ((fnow > enow) == (fnow > -enow))
    {
      g0 = enow;
      eindex++;
      if (eindex < e.size())
	{
	  enow = e[eindex];
	}
    }
  else
    {
      g0 = fnow;
      findex++;
      if (findex < f.size())
	{
	  fnow = f[findex];
	}
    }
  if ((eindex <e.size()) && ((findex >= f.size()) || ((fnow > enow) == (fnow > -enow))))
    {
      fastTwoSum(enow, g0, Qnew, q);
      eindex++;
      if (eindex < e.size())
	{
	  enow = e[eindex];
	}
    }
  else
    {
      fastTwoSum(fnow, g0, Qnew, q);
      findex++;
      if (findex <f.size())
	{
	  fnow = f[findex];
	}
    }
  Q = Qnew;
  for (size_t count = 2, endp=e.size() + f.size(); count < endp; count++)
    {
      if ((eindex < e.size()) && ((findex >= f.size()) || ((fnow > enow) == (fnow > -enow))))
	{
	  fastTwoSum(enow, q, R, hh);
	  eindex++;
	  if (eindex < e.size())
	    {
	      enow = e[eindex];
	    }
	}
      else
	{
	  fastTwoSum(fnow, q, R, hh);
	  findex++;
	  if (findex < f.size())
	    {
	      fnow = f[findex];
	    }
	}
      twoSum(Q, R, Qnew, q);
      Q = Qnew;
      if (fabs(hh) > 0)
	{
	  ans.push_back(hh);
	}
    }
  if (fabs(q) > 0) {
    ans.push_back(q);
  }
  if ((Q != 0.0) || ans.empty()) {
    ans.push_back(Q);
  }
  return ans;
}

void scaleExpansionZeroElim(vector<double> const& e, double b,vector<double>
			    &ans)
{
  ans.clear();
  ans.reserve(10);
  double Q, hh;
  double product0, product1,sum;
  int eindex;

  twoProduct(e[0], b, Q, hh);
  if (fabs(hh) > 0)
    {
      ans.push_back(hh);
    }
  int n=static_cast<int>(e.size());
  for (eindex = 1; eindex<n; ++eindex)
    {
      twoProduct(e[static_cast<size_t>(eindex)], b, product1, product0);
      twoSum(Q, product0,sum, hh);
      if (fabs(hh) > 0)
	{
	  ans.push_back(hh);
	}
      fastTwoSum(product1,sum, Q, hh);
      if (fabs(hh) > 0)
	{
	  ans.push_back(hh);
	}
    }
  if ((Q != 0.0) || ans.empty()) {
    ans.push_back(Q);
  }
}

vector<double> compress(vector<double> const& e)
{
  vector<double> temp, ans;
  ans.reserve(8);

  int bottom = static_cast<int>(e.size()) - 1;
  double Q = e[static_cast<size_t>(bottom)];
  for (int eindex =static_cast<int>(e.size()) - 2; eindex >= 0; eindex--)
    {
      const double enow = e[static_cast<size_t>(eindex)];
      double q, Qnew;
      fastTwoSum(Q, enow, Qnew, q);
      if (fabs(q) > 0)
	{
	  temp.insert(temp.begin(), Qnew);
	  bottom--;
	  Q = q;
	}
      else
	{
	  Q = Qnew;
	}
    }
  for (size_t i = static_cast<size_t>(bottom) + 1; i <e.size(); i++)
    {
      const double current = temp[i];
      double q, Qnew;
      fastTwoSum(current, Q, Qnew, q);
      if (fabs(q) > 0)
	{
	  ans.push_back(q);
	}
      Q = Qnew;
    }
  ans.push_back(Q);
  return ans;
}

double estimate(vector<double> const& e)
{
  double Q;

  Q = e[0];
  for (size_t i = 1; i < e.size(); i++)
    {
      Q += e[i];
    }
  return Q;
}
