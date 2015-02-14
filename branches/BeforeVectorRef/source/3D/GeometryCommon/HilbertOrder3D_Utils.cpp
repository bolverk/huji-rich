#include <cmath>
#include "HilbertOrder3D_Utils.hpp"

int EstimateHilbertIterationNum(vector<Vector3D> const& cor)
{
	return (int)ceil(log(pow(static_cast<double>(cor.size()), (1.0 / 3.0))) / log(2.0));
}


// Adjust the vector of points so all coordinate values will be in the range [0,1]
void AdjustPoints(vector<Vector3D> const & vPointsIn, vector<Vector3D> & vPointsOut)
{
	// The output vector:
	vPointsOut.resize(vPointsIn.size());
	// Vectors holding the X,Y,Z coordinates:
	vector<double> vPointsX, vPointsY, vPointsZ;
	
	Split(vPointsIn, vPointsX, vPointsY, vPointsZ);

	// TOOD - what if all points have the same x (or y, or z) values?

	double dbMinX = *min_element(vPointsX.begin(), vPointsX.end());
	double dbMinY = *min_element(vPointsY.begin(), vPointsY.end());
	double dbMinZ = *min_element(vPointsZ.begin(), vPointsZ.end());

	double dbMaxX = *max_element(vPointsX.begin(), vPointsX.end());
	double dbMaxY = *max_element(vPointsY.begin(), vPointsY.end());
	double dbMaxZ = *max_element(vPointsZ.begin(), vPointsZ.end());
	// The scale factor:
	double dbScaleX = dbMaxX - dbMinX;
	double dbScaleY = dbMaxY - dbMinY;
	double dbScaleZ = dbMaxZ - dbMinZ;
	// To prevent division by zero (very unlikely - double precision!)
	bool bFlagX = dbScaleX == 0;
	bool bFlagY = dbScaleY == 0;
	bool bFlagZ = dbScaleZ == 0;

	// X coordinate:
	if (!bFlagX)
	{
		for (size_t ii = 0; ii < vPointsIn.size(); ++ii)
		{
			// Scale the X coordinate:
			vPointsX[ii] -= dbMinX;
			vPointsX[ii] /= dbScaleX;
		}
	}
	else
	{
		fill(vPointsX.begin(), vPointsX.end(), 0);
	}
	// Y coordinate:
	if (!bFlagY)
	{
		for (size_t ii = 0; ii < vPointsIn.size(); ++ii)
		{
			// Scale the X coordinate:
			vPointsY[ii] -= dbMinY;
			vPointsY[ii] /= dbScaleY;
		}
	}
	else
	{
		fill(vPointsY.begin(), vPointsY.end(), 0);
	}
	// Z coordinate:
	if (!bFlagZ)
	{
		for (size_t ii = 0; ii < vPointsIn.size(); ++ii)
		{
			vPointsZ[ii] -= dbMinZ;
			vPointsZ[ii] /= dbScaleZ;
		}
	}
	else
	{
		fill(vPointsZ.begin(), vPointsZ.end(), 0);
	}

	// Store back the coordiantes in a single output vector:
	for (size_t ii = 0; ii < vPointsIn.size(); ++ii)
	{
		vPointsOut[ii].x = vPointsX[ii];
		vPointsOut[ii].y = vPointsY[ii];
		vPointsOut[ii].z = vPointsZ[ii];
	}

	return;
}

void FindEqualIndices(vector<unsigned long long int> const & vD_sorted, vector<vector<size_t> > & vOut)
{
	vector<unsigned long long int> vD_sorted_cpy = vD_sorted;
	vector<unsigned long long int> vD_sorted_unq = vD_sorted;
	//vector<unsigned long long int>::iterator it = adjacent_find(vD_sorted.begin(), vD_sorted.end());

	vector<unsigned long long int>::iterator it1, itPrev, itCur;
	it1 = unique(vD_sorted_unq.begin(), vD_sorted_unq.end());

	vD_sorted_unq.resize(distance(vD_sorted_unq.begin(), it1));
	
	if (vD_sorted.size() == vD_sorted_unq.size())
	{
		return;
	}

	vOut.reserve(vD_sorted.size() - vD_sorted_unq.size());

	int iCurPrevDist = 0;

	itPrev = vD_sorted_cpy.begin();
	for (it1 = vD_sorted_unq.begin()+1; it1 != vD_sorted_unq.end(); ++it1)
	{
		if (distance( it1 , vD_sorted_unq.end() ) == 0)
		{
			itCur = vD_sorted_cpy.end();
		}
		else
		{
			itCur = find(itPrev, vD_sorted_cpy.end(), *it1);
		}

		iCurPrevDist = static_cast<int>(distance(itPrev, itCur));
		if (1 < iCurPrevDist)
		{
		  int iBase = static_cast<int>(distance(vD_sorted_cpy.begin(), itPrev));
			vector<size_t> vInd( iCurPrevDist );
			// C++11
			// iota(vInd.begin(), vInd.end(), iBase);
			for (int ii = 0; ii < iCurPrevDist; ++ii)
			{
				vInd[ii] = iBase + ii;
			}
			vOut.push_back(vInd);
		}

		itPrev = itCur;
	}

	iCurPrevDist = static_cast<int>(distance(itPrev, vD_sorted_cpy.end()));
	if (1 < iCurPrevDist )
	{
	  int iBase = static_cast<int>(distance(vD_sorted_cpy.begin(), itPrev));

		vector<size_t> vInd(iCurPrevDist);
		for (int ii = 0; ii < iCurPrevDist; ++ii)
		{
			vInd[ii] = iBase + ii;
		}
		vOut.push_back(vInd);
	}

	//	int x = 0;

	return;
}
