#include "HilbertOrder3D.hpp"
#include "HilbertOrder3D_Utils.hpp"
#include <vector>
#include <algorithm>
#include <ctime>
#include <iostream>

#define NUMBER_OF_SHAPES 24
#define MAX_ROTATION_LENGTH 5
#define PI 3.14159

//! \brief The elementary Hilbert Curve shape
class HilbertCurve3D_shape
{
public:
	//! \brief Class constructor
	HilbertCurve3D_shape();
	//! \brief Comparison
  //! \param shape of hilbert curve
  //! \return True if two shapes are the same
	bool operator==(HilbertCurve3D_shape & shape);
	//! \brief An array of the 7 unit vector steps defining the shape
	vector<Vector3D> m_vShapePoints;

};
// Constructor - uses the reference Hilbert curve shape:
HilbertCurve3D_shape::HilbertCurve3D_shape():
  m_vShapePoints(vector<Vector3D> ())
{
	m_vShapePoints.resize(7);
	m_vShapePoints[0] = Vector3D(0, 0, -1);
	m_vShapePoints[1] = Vector3D(0, 1, 0);
	m_vShapePoints[2] = Vector3D(0, 0, 1);
	m_vShapePoints[3] = Vector3D(-1, 0, 0);
	m_vShapePoints[4] = Vector3D(0, 0, -1);
	m_vShapePoints[5] = Vector3D(0, -1, 0);
	m_vShapePoints[6] = Vector3D(0, 0, 1);
}

// Compare to a given Hilbert curve shape by comparing pairs of shape points:
bool HilbertCurve3D_shape::operator==(HilbertCurve3D_shape & shape)
{
	bool b = true;
	for (size_t ii = 0; ii < m_vShapePoints.size(); ++ii)
	{
		b = b && (m_vShapePoints[ii] == shape.m_vShapePoints[ii]);
	}

	return b;
}

//! \brief Hilbert Curve
class HilbertCurve3D
{
public:
	//! \brief Constructor
	HilbertCurve3D(void);
  /*! \brief Calculate the Hilbert curve distance of a given point, given a required number of iterations
    \param rvPoint Pivot
    \param numOfIterations Number of iterations
    \return Position index
   */
	size_t Hilbert3D_xyz2d(Vector3D const & rvPoint, int numOfIterations);

private:
	// Rotate a shape according to a given rotation scheme (in-place):
	void RotateShape(int iShapeIndex, vector<int> vAxes);
	// Rotate a shape according to rotation index, and return the rotated shape:
	void RotateShape(HilbertCurve3D_shape const & roShape, HilbertCurve3D_shape & roShapeOut, int iRotationIndex);
	/*!
	\brief Returns the rotation scheme, according to a rotation index
	\param piRotation - a pointer to the output rotation scheme vector
	\param iRotationIndex - the desired rotation index
	\return The rotation scheme length (the size of the array given by piRotation)
	*/
	int GetRotation(int * piRotation, int iRotationIndex);
	// Find the index of a given shape object:
	int FindShapeIndex(HilbertCurve3D_shape & roShape);
	// Create the recursion rule:
	void BuildRecursionRule();
	// Create the shape order, for all shapes (the order of octants):
	void BuildShapeOrder();

	// Stores all rotated shapes:
	vector<HilbertCurve3D_shape> m_vRotatedShapes;
	// Stores all rotation schemes:
	vector < vector<int> > m_vRotations;

	// An array of the 8 integers defining the recursion rule of the shape:
	vector< vector<int> > m_vShapeRecursion;

	// A 2x2x2 matrix indicating the 3 dimensional shape order
	// array< array<int , 8 > , NUMBER_OF_SHAPES > m_mShapeOrder;
	int m_mShapeOrder[NUMBER_OF_SHAPES][2][2][2];
};

// Constructor - performs all required initiallizations and preprocessing:
HilbertCurve3D::HilbertCurve3D():
  m_vRotatedShapes(),
  m_vRotations(),
  m_vShapeRecursion()
{
	m_vRotatedShapes.resize(NUMBER_OF_SHAPES);
	m_vRotations.resize(NUMBER_OF_SHAPES);
	m_vShapeRecursion.resize(NUMBER_OF_SHAPES);
	for(size_t i=0;i<NUMBER_OF_SHAPES;++i)
		m_vShapeRecursion[i].resize(8);
	int rot[MAX_ROTATION_LENGTH];
	for (int iRotIndex = 1; iRotIndex < NUMBER_OF_SHAPES; ++iRotIndex)
	{
	  const int iRotLength = GetRotation(rot, iRotIndex);
	  m_vRotations[static_cast<size_t>(iRotIndex)].assign(rot, rot + iRotLength);
	}

	for (int ii = 1; ii < NUMBER_OF_SHAPES; ++ii)
	{
	  RotateShape(ii, m_vRotations[static_cast<size_t>(ii)]);
	}

	BuildRecursionRule();
	BuildShapeOrder();
}

// FindShapeIndex - returns the index of a shape:
int HilbertCurve3D::FindShapeIndex(HilbertCurve3D_shape & roShape)
{
	for (int ii = 0; ii < NUMBER_OF_SHAPES; ++ii)
	{
	  if (roShape == m_vRotatedShapes[static_cast<size_t>(ii)])
		{
			return ii;
		}
	}
	// TODO - manage this kind of return value (error)
	return -1;
}

void HilbertCurve3D::BuildRecursionRule()
{
	// Reference recursion rule:
	m_vShapeRecursion[0][0] = 12;
	m_vShapeRecursion[0][1] = 16;
	m_vShapeRecursion[0][2] = 16;
	m_vShapeRecursion[0][3] = 2;
	m_vShapeRecursion[0][4] = 2;
	m_vShapeRecursion[0][5] = 14;
	m_vShapeRecursion[0][6] = 14;
	m_vShapeRecursion[0][7] = 10;

	HilbertCurve3D_shape oTempShape;
	// What about ii=0? not necessary 
	for (int ii = 0; ii < NUMBER_OF_SHAPES; ++ii)
	{
		for (int jj = 0; jj < 8; ++jj)
		{
			// Rotate the appropriate block of the reference recursion rule, according to the ii rotation scheme:
		  RotateShape(m_vRotatedShapes[static_cast<size_t>(m_vShapeRecursion[0][static_cast<size_t>(jj)])], oTempShape, ii);
			// Find the shape index of the rotated shape:
		  m_vShapeRecursion[static_cast<size_t>(ii)][static_cast<size_t>(jj)] = FindShapeIndex(oTempShape);
		}
	}

	return;
}

// Return the rotation scheme of the ii rotation (manually precalculated):
int HilbertCurve3D::GetRotation(int * piRotation, int iRotationIndex)
{
	switch (iRotationIndex)
	{
	case 1:
		piRotation[0] = 1;
		return 1;
	case 2:
		piRotation[0] = 1;
		piRotation[1] = 1;
		return 2;
	case 3:
		piRotation[0] = 1;
		piRotation[1] = 1;
		piRotation[2] = 1;
		return 3;
	
	case 4:
		piRotation[0] = 2;
		return 1;
	case 5:
		piRotation[0] = 2;
		piRotation[1] = 2;
		return 2;
	case 6:
		piRotation[0] = 2;
		piRotation[1] = 2;
		piRotation[2] = 2;
		return 3;

	case 7:
		piRotation[0] = 3;
		return 1;
	case 8:
		piRotation[0] = 3;
		piRotation[1] = 3;
		return 2;
	case 9:
		piRotation[0] = 3;
		piRotation[1] = 3;
		piRotation[2] = 3;
		return 3;

	case 10:
		piRotation[0] = 1;
		piRotation[1] = 2;
		return 2;
	case 11:
		piRotation[0] = 2;
		piRotation[1] = 1;
		return 2;
	case 12:
		piRotation[0] = 3;
		piRotation[1] = 1;
		return 2;
	case 13:
		piRotation[0] = 2;
		piRotation[1] = 3;
		return 2;

	case 14:
		piRotation[0] = 1;
		piRotation[1] = 1;
		piRotation[2] = 1;
		piRotation[3] = 3;
		return 4;
	case 15:
		piRotation[0] = 2;
		piRotation[1] = 2;
		piRotation[2] = 2;
		piRotation[3] = 1;
		return 4;
	case 16:
		piRotation[0] = 3;
		piRotation[1] = 3;
		piRotation[2] = 3;
		piRotation[3] = 2;
		return 4;
	case 17:
		piRotation[0] = 3;
		piRotation[1] = 2;
		piRotation[2] = 3;
		piRotation[3] = 2;
		return 4;
	case 18:
		piRotation[0] = 3;
		piRotation[1] = 2;
		piRotation[2] = 2;
		return 3;

	case 19:
		piRotation[0] = 1;
		piRotation[1] = 1;
		piRotation[2] = 1;
		piRotation[3] = 2;
		piRotation[4] = 2;
		return 5;
	case 20:
		piRotation[0] = 3;
		piRotation[1] = 3;
		piRotation[2] = 2;
		return 3;
	case 21:
		piRotation[0] = 3;
		piRotation[1] = 3;
		piRotation[2] = 1;
		return 3;
	case 22:
		piRotation[0] = 2;
		piRotation[1] = 2;
		piRotation[2] = 3;
		return 3;
	case 23:
		piRotation[0] = 1;
		piRotation[1] = 1;
		piRotation[2] = 2;
		return 3;

	default:
		return 0;
		//		break;
	}
}

// Rotate a shape:
void HilbertCurve3D::RotateShape(int iShapeIndex, vector<int> vAxes)
{
	int iSign = 0;

	for (size_t ii = 0; ii < 7; ++ii)
	{
		for (size_t iAx = 0; iAx < vAxes.size(); ++iAx)
		{
			// A trick to find the sign of vAxes[iAx]:
			iSign = (vAxes[iAx] > 0) - (vAxes[iAx] < 0);

			switch (abs(vAxes[iAx]))
			{
			case 1:
			  m_vRotatedShapes[static_cast<size_t>(iShapeIndex)].m_vShapePoints[ii].RotateX( iSign * PI / 2 );
				break;
			case 2:
			  m_vRotatedShapes[static_cast<size_t>(iShapeIndex)].m_vShapePoints[ii].RotateY( iSign * PI / 2 );
				break;
			case 3:
			  m_vRotatedShapes[static_cast<size_t>(iShapeIndex)].m_vShapePoints[ii].RotateZ( iSign * PI / 2 );
				break;
			default:
				break;
			}
		}
		// Round off the results:
		m_vRotatedShapes[static_cast<size_t>(iShapeIndex)].m_vShapePoints[ii].Round();
	}
}

void HilbertCurve3D::RotateShape(HilbertCurve3D_shape const & roShape, HilbertCurve3D_shape & roShapeOut , int iRotationIndex)
{
	int iSign = 0;

	vector<int> vAxes = m_vRotations[static_cast<size_t>(iRotationIndex)];
	roShapeOut = roShape;

	for (int ii = 0; ii < 7; ++ii)
	{
		for (size_t iAx = 0; iAx < vAxes.size(); ++iAx)
		{
			// A trick to find the sign of vAxes[iAx]:
			iSign = (vAxes[iAx] > 0) - (vAxes[iAx] < 0);

			switch (abs(vAxes[iAx]))
			{
			case 1:
			  roShapeOut.m_vShapePoints[static_cast<size_t>(ii)].RotateX(iSign * PI / 2);
				break;
			case 2:
			  roShapeOut.m_vShapePoints[static_cast<size_t>(ii)].RotateY(iSign * PI / 2);
				break;
			case 3:
			  roShapeOut.m_vShapePoints[static_cast<size_t>(ii)].RotateZ(iSign * PI / 2);
				break;
			default:
				break;
			}
		}
		// Round off the results:
		roShapeOut.m_vShapePoints[static_cast<size_t>(ii)].Round();
	}

	return;
}

void HilbertCurve3D::BuildShapeOrder()
{
	vector<int> vShapeVerticesX(8);
	vector<int> vShapeVerticesY(8);
	vector<int> vShapeVerticesZ(8);

	vShapeVerticesX[0] = 0;
	vShapeVerticesY[0] = 0;
	vShapeVerticesZ[0] = 0;

	for (size_t iShapeInd = 0; iShapeInd < NUMBER_OF_SHAPES; ++iShapeInd)
	{
		for (size_t ii = 0; ii < m_vRotatedShapes[iShapeInd].m_vShapePoints.size(); ++ii)
		{
			vShapeVerticesX[ii + 1] = static_cast<int>(vShapeVerticesX[ii] + m_vRotatedShapes[iShapeInd].m_vShapePoints[ii].x);
			vShapeVerticesY[ii + 1] = static_cast<int>(vShapeVerticesY[ii] + m_vRotatedShapes[iShapeInd].m_vShapePoints[ii].y);
			vShapeVerticesZ[ii + 1] = static_cast<int>(vShapeVerticesZ[ii] + m_vRotatedShapes[iShapeInd].m_vShapePoints[ii].z);
		}

		int iMinX = *min_element(vShapeVerticesX.begin(), vShapeVerticesX.end());
		int iMinY = *min_element(vShapeVerticesY.begin(), vShapeVerticesY.end());
		int iMinZ = *min_element(vShapeVerticesZ.begin(), vShapeVerticesZ.end());

		for (size_t jj = 0; jj < vShapeVerticesX.size(); ++jj)
		{
			vShapeVerticesX[jj] -= iMinX;
			vShapeVerticesY[jj] -= iMinY;
			vShapeVerticesZ[jj] -= iMinZ;
		}

		for (size_t kk = 0; kk < vShapeVerticesX.size(); ++kk)
		{
		  m_mShapeOrder[iShapeInd][vShapeVerticesX[kk]][vShapeVerticesY[kk]][vShapeVerticesZ[kk]] = static_cast<int>(kk);
		}
	}

	return;
}

size_t HilbertCurve3D::Hilbert3D_xyz2d(Vector3D const & rvPoint, int numOfIterations)
{
	// Extract the coordinates:
	double x = rvPoint.x;
	double y = rvPoint.y;
	double z = rvPoint.z;

	// The output distance along the 3D-Hilbert Curve:
	size_t d = 0;

	// The current shape index:
	int iCurrentShape = 0;
	for (int iN = 1; iN <= numOfIterations; ++iN)
	{
		// Calculate the current power of 0.5:
	  const double dbPow2 = (static_cast<double>(1)) / (1 << iN);
	  const bool bX = x > dbPow2;
	  const bool bY = y > dbPow2;
	  const bool bZ = z > dbPow2;

		x -= dbPow2*bX;
		y -= dbPow2*bY;
		z -= dbPow2*bZ;

		// Multiply the distance by 8 (for every recursion iteration):
		d = d << 3;
	const int iOctantNum = m_mShapeOrder[iCurrentShape][bX][bY][bZ];
	d = d + static_cast<size_t>(iOctantNum);
	iCurrentShape = m_vShapeRecursion[static_cast<size_t>(iCurrentShape)][static_cast<size_t>(iOctantNum)];
	}

	//	int a = 0;
	return d;
}

vector<size_t> HilbertOrder3D(vector<Vector3D> const& cor)
{
	// If only 1 or 2 points are provided - do not reorder them
	if ( 2 >= cor.size() )
	{
		vector<size_t> vIndSort( cor.size() );
		for (size_t ii = 0; ii < cor.size(); ++ii)
		{
			vIndSort[ii] = ii;
		}
		return vIndSort;
	}
	// Create a 3D-Hilbert Curve Object:
	HilbertCurve3D oHilbert;

	// Allocate an output vector:
	int N = static_cast<int>(cor.size());
	vector<size_t> vOut;
	vOut.reserve(cor.size());
	
	// Estimate the number of required iterations:
	int numOfIterations = EstimateHilbertIterationNum(cor);

	vector<Vector3D> vAdjustedPoints;

	// Adjust the points coordinates to the unit cube:
	AdjustPoints(cor, vAdjustedPoints);

	// Run throught the points, and calculate the Hilbert distance of each:
	for (int ii = 0; ii < N; ++ii)
	{
	  vOut.push_back(oHilbert.Hilbert3D_xyz2d(vAdjustedPoints[static_cast<size_t>(ii)], numOfIterations+6));
		//vOut.push_back(oHilbert.Hilbert3D_xyz2d(vAdjustedPoints[ii], 2));
	}
	// Get the sorting indices:
	vector<size_t> vIndSort;
	sort_index(vOut,vIndSort);
	// Reorder the Hilbert distances vector (according to the sorting indices):
	reorder( vOut, vIndSort );

	// Find indices with repeated Hilbert distance:
	vector<vector<size_t> > vEqualIndices;
	FindEqualIndices(vOut, vEqualIndices);
	
	// If all points have different Hilbert distances, return the sorting indices:
	if (vEqualIndices.empty())
	{
		return vIndSort;
	}
	else
	{
		for (size_t ii = 0; ii < vEqualIndices.size(); ++ii)
		{
			vector<Vector3D> vPointsInner(vEqualIndices[ii].size() );
			vector<size_t> vIndInner(vEqualIndices[ii].size());
			vector<size_t> vIndSortInner(vEqualIndices[ii].size());
			vector<size_t> vIndSortInner_cpy(vEqualIndices[ii].size());

			// Store the points with the equal indices
			for (size_t jj = 0; jj < vEqualIndices[ii].size() ; ++jj)
			{
				vIndInner[jj] = vIndSort[vEqualIndices[ii][jj]];
				vPointsInner[jj] = cor[vIndInner[jj]];
				vIndSortInner_cpy[jj] = vIndSort[vIndInner[jj]];
			}
			
			// Sort the repeated points:
			vIndSortInner = HilbertOrder3D(vPointsInner);
	//		vector<size_t> vIndSortTemp = vIndSort;
			for (size_t kk = 0; kk < vIndSortInner.size(); ++kk)
			{
				vIndSort[vIndInner[kk]] = vIndSortInner_cpy[vIndSortInner[kk]];
			}
		}

		// Return the sorting indices:
		return vIndSort;
	}
}
