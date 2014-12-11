#include "geotests.hpp"
double const epsilon = 1.1102230246251565e-016;

double orient2dAdapt(const TripleConstRef<Vector2D>& points, double detsum)
{
	boost::array<double,2> A1,A2;
	vector<double> B1, B2;
	vector<double> C;
	double acx, acy, bcx, bcy, acxtail, acytail, bcxtail, bcytail;
	double detLeft, detRight, detLeftErr, detRightErr, det, err;
	double s1, t1, s0, t0;
	// Error ranges for different accuracies.
	double resultErrBound = (3.0 + 8.0 * epsilon) * epsilon;
	double errBoundB = (2.0 + 12.0 * epsilon) * epsilon;
	double errBoundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;

	twoDiff(points.first.x, points.third.x, acx, acxtail);
	twoDiff(points.second.x, points.third.x, bcx, bcxtail);
	twoDiff(points.first.y, points.third.y, acy, acytail);
	twoDiff(points.second.y, points.third.y, bcy, bcytail);

	twoProduct(acx, bcy, detLeft, detLeftErr);
	twoProduct(acy, bcx, detRight, detRightErr);
	A1[0] = detLeft;
	A1[1] = detLeftErr;
	A2[0] = detRight;
	A2[1] = detRightErr;
	B1 = twoTwoDiff(A1, A2);
	det = estimate(B1);
	err = errBoundB * detsum;
	// Several cases in which precision is enough.
	if ((det >= err) || (-det >= err))
	{
		return det;
	}
	if ((acxtail == 0.0) && (acytail == 0.0) && (bcxtail == 0.0) && (bcytail == 0.0))
	{
		return det;
	}
	err = errBoundC * detsum + resultErrBound * fabs(det);
	det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
	if ((det >= err) || (-det >= err))
	{
		return det;
	}

	twoProduct(acxtail, bcy, s1, s0);
	twoProduct(acytail, bcx, t1, t0);
	A1[0] = s1;
	A1[1] = s0;
	A2[0] = t1;
	A2[1] = t0;
	B2 = twoTwoDiff(A1, A2);
	C = fastExpansionSumZeroElim(B1, B2);

	twoProduct(acx, bcytail, s1, s0);
	twoProduct(acy, bcxtail, t1, t0);
	A1[0] = s1;
	A1[1] = s0;
	A2[0] = t1;
	A2[1] = t0;
	B1 = twoTwoDiff(A1, A2);
	C = fastExpansionSumZeroElim(C, B1);

	twoProduct(acxtail, bcytail, s1, s0);
	twoProduct(acytail, bcxtail, t1, t0);
	A1[0] = s1;
	A1[1] = s0;
	A2[0] = t1;
	A2[1] = t0;
	B1 = twoTwoDiff(A1, A2);
	C = fastExpansionSumZeroElim(C, B1);
	return(C.back());
}

double orient2d(const TripleConstRef<Vector2D>& points)
{
	double detleft, detright, det;
	double detsum, errbound;
	double errBoundA = (3.0 + 16.0 * epsilon) * epsilon; // Acceptable error range for this calculation precision.

	// Calculating a 2x2 determinant.
	detleft = (points.first.x - points.third.x) * (points.second.y - points.third.y);
	detright = (points.first.y - points.third.y) * (points.second.x - points.third.x);
	det = detleft - detright;

	if (detleft > 0.0)
	{
		if (detright <= 0.0) // Determinant sign is certain.
		{
			return det;
		}
		else
		{
			detsum = detleft + detright;
		}
	}
	else
	{
		if (detleft < 0.0)
		{
			if (detright >= 0.0) //Determinant sign is certain.
			{
				return det;
			}
			else
			{
				detsum = -detleft - detright;
			}
		}
		else
		{
			return det;
		}
	}
	errbound = errBoundA * detsum;
	if ((det >= errbound) || (-det >= errbound)) // Calculation in in acceptable error range.
	{
		return det;
	}
	return orient2dAdapt(points, detsum);
}

double incircleadapt(const Vector2D& point_1,
		     const Vector2D& point_2,
		     const Vector2D& point_3,
		     const Vector2D& point_4,
		     double permanent)
{
	double adx, bdx, cdx, ady, bdy, cdy;
	double adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
	double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
	double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
	double res1, res2, err1, err2;
	double det, errbound;
	// Several error ranges for different calculation precisions.
	double errBoundB = (4.0 + 48.0 * epsilon) * epsilon;
	double errBoundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
	double resultErrBound = (3.0 + 8.0 * epsilon) * epsilon;

	boost::array<double,2> A1, A2;
	vector<double> B1, B2, B3, C1, C2, C3, D1, D2, D3, E1, E2, E3, F;
	vector<double> tmp1, tmp2;

	twoDiff(point_1.x, point_4.x, adx, adxtail);
	twoDiff(point_2.x, point_4.x, bdx, bdxtail);
	twoDiff(point_3.x, point_4.x, cdx, cdxtail);
	twoDiff(point_1.y, point_4.y, ady, adytail);
	twoDiff(point_2.y, point_4.y, bdy, bdytail);
	twoDiff(point_3.y, point_4.y, cdy, cdytail);

	twoProduct(bdx, cdy, bdxcdy1, bdxcdy0);
	twoProduct(cdx, bdy, cdxbdy1, cdxbdy0);
	A1[0] = bdxcdy1;
	A1[1] = bdxcdy0;
	A2[0] = cdxbdy1;
	A2[1] = cdxbdy0;
	B1 = twoTwoDiff(A1, A2);
	scaleExpansionZeroElim(B1,adx,C1);
	scaleExpansionZeroElim(C1,adx,tmp1);
	scaleExpansionZeroElim(B1, ady,C1);
	scaleExpansionZeroElim(C1,ady,tmp2);
	C1 = fastExpansionSumZeroElim(tmp1,tmp2);

	twoProduct(cdx, ady, cdxady1, cdxady0);
	twoProduct(adx, cdy, adxcdy1, adxcdy0);
	A1[0] = cdxady1;
	A1[1] = cdxady0;
	A2[0] = adxcdy1;
	A2[1] = adxcdy0;
	B2 = twoTwoDiff(A1, A2);
	scaleExpansionZeroElim(B2,bdx,tmp2);
	scaleExpansionZeroElim(tmp2,bdx,tmp1);
	scaleExpansionZeroElim(B2,bdy,C2);
	scaleExpansionZeroElim(C2,bdy,tmp2);
	C2 = fastExpansionSumZeroElim(tmp1,tmp2);

	twoProduct(adx, bdy, adxbdy1, adxbdy0);
	twoProduct(bdx, ady, bdxady1, bdxady0);
	A1[0] = adxbdy1;
	A1[1] = adxbdy0;
	A2[0] = bdxady1;
	A2[1] = bdxady0;
	B3 = twoTwoDiff(A1, A2);
	scaleExpansionZeroElim(B3,cdx,tmp2);
	scaleExpansionZeroElim(tmp2,cdx,tmp1);
	scaleExpansionZeroElim(B3,cdy,C3);
	scaleExpansionZeroElim(C3,cdy,tmp2);
	C3 = fastExpansionSumZeroElim(tmp1, tmp2);

	F = fastExpansionSumZeroElim(C1, C2);
	F = fastExpansionSumZeroElim(F, C3);
	det = estimate(F); //Estimating determinant.
	// Several cases in which precision is enough.
	errbound = errBoundB * permanent;
	if ((det >= errbound) || (-det >= errbound))
	{
		return det;
	}
	if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
		&& (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0))
	{
			return det;
	}
	errbound = errBoundC * permanent + resultErrBound * fabs(det);
	det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail)
		- (bdy * cdxtail + cdx * bdytail)) + 2.0 * (adx * adxtail + ady * adytail)
		* (bdx * cdy - bdy * cdx)) + ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail)
		- (cdy * adxtail + adx * cdytail)) + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
		+ ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
		+ 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
	if ((det >= errbound) || (-det >= errbound))
	{
		return det;
	}

	// A more precise calculation.
	if ((bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) || (cdytail != 0.0))
	{
	  double adxadx1;
	  double adyady1;
	  double adxadx0;
	  double adyady0;
		square(adx, adxadx1, adxadx0);
		square(ady, adyady1, adyady0);
		A1[0] = adxadx1;
		A1[1] = adxadx0;
		A2[0] = adyady1;
		A2[1] = adyady0;
		C1 = twoTwoSum(A1, A2);
	}
	if ((cdxtail != 0.0) || (cdytail != 0.0) || (adxtail != 0.0) || (adytail != 0.0))
	{
	  double bdxbdx1;
	  double bdybdy1;
	  double bdxbdx0;
	  double bdybdy0;
		square(bdx, bdxbdx1, bdxbdx0);
		square(bdy, bdybdy1, bdybdy0);
		A1[0] = bdxbdx1;
		A1[1] = bdxbdx0;
		A2[0] = bdybdy1;
		A2[1] = bdybdy0;
		C2 = twoTwoSum(A1, A2);
	}
	if ((adxtail != 0.0) || (adytail != 0.0) || (bdxtail != 0.0) || (bdytail != 0.0))
	{
	  double cdxcdx1;
	  double cdycdy1;
	  double cdxcdx0;
	  double cdycdy0;
		square(cdx, cdxcdx1, cdxcdx0);
		square(cdy, cdycdy1, cdycdy0);
		A1[0] = cdxcdx1;
		A1[1] = cdxcdx0;
		A2[0] = cdycdy1;
		A2[1] = cdycdy0;
		C3 = twoTwoSum(A1, A2);
	}
	if (adxtail != 0.0)
	{
		vector<double> helper;
		scaleExpansionZeroElim(B1,adxtail,D1);
		scaleExpansionZeroElim(D1,2.0*adx,tmp1);
		scaleExpansionZeroElim(C3, adxtail,helper);
		scaleExpansionZeroElim(helper,bdy,tmp2);
		tmp2 = fastExpansionSumZeroElim(tmp1,tmp2);
		scaleExpansionZeroElim(C2, adxtail,helper);
		scaleExpansionZeroElim(helper,-cdy,tmp1);
		tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
		F = fastExpansionSumZeroElim(F, tmp1);
	}
	if (adytail != 0.0)
	{
		scaleExpansionZeroElim(B1, adytail,D2);
		scaleExpansionZeroElim(D2,2.0*ady,tmp1);
		vector<double> helper;
		scaleExpansionZeroElim(C2,adytail,tmp2);
		scaleExpansionZeroElim(tmp2,cdx,helper);
		tmp2 = fastExpansionSumZeroElim(tmp1,helper);
		scaleExpansionZeroElim(C3,adytail,tmp1);
		scaleExpansionZeroElim(tmp1,-bdx,helper);
		tmp1 = fastExpansionSumZeroElim(helper,tmp2);
		F = fastExpansionSumZeroElim(F,tmp1);
	}
	if (bdxtail != 0.0)
	{
		vector<double> helper;
		scaleExpansionZeroElim(B2, bdxtail,D3);
		scaleExpansionZeroElim(D3, 2.0 * bdx,tmp1);
		scaleExpansionZeroElim(C1, bdxtail,tmp2);
		scaleExpansionZeroElim(tmp2, cdy,helper);
		tmp2 = fastExpansionSumZeroElim(tmp1,helper);
		scaleExpansionZeroElim(C3, bdxtail,tmp1);
		scaleExpansionZeroElim(tmp1, -ady,helper);
		tmp1 = fastExpansionSumZeroElim(helper, tmp2);
		F = fastExpansionSumZeroElim(F, tmp1);
	}
	if (bdytail != 0.0)
	{
		vector<double> helper;
		scaleExpansionZeroElim(B2, bdytail,B1);
		scaleExpansionZeroElim(B1, 2.0 * bdy,tmp1);
		scaleExpansionZeroElim(C3, bdytail,helper);
		scaleExpansionZeroElim(helper, adx,tmp2);
		tmp2 = fastExpansionSumZeroElim(tmp1, tmp2);
		scaleExpansionZeroElim(C1, bdytail,helper);
		scaleExpansionZeroElim(helper,-cdx,tmp1);
		tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
		F = fastExpansionSumZeroElim(F, tmp1);
	}
	if (cdxtail != 0.0)
	{
		vector<double> helper;
		scaleExpansionZeroElim(B3, cdxtail,B2);
		scaleExpansionZeroElim(B2, 2.0 * cdx,tmp1);
		scaleExpansionZeroElim(C2, cdxtail,helper);
		scaleExpansionZeroElim(helper, ady,tmp2);
		tmp2 = fastExpansionSumZeroElim(tmp1, tmp2);
		scaleExpansionZeroElim(C1, cdxtail,helper);
		scaleExpansionZeroElim(helper, -bdy,tmp1);
		tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
		F = fastExpansionSumZeroElim(F, tmp1);
	}
	if (cdytail != 0.0)
	{
		vector<double> helper;
		scaleExpansionZeroElim(B3,cdytail,helper);
		B3=helper;
		scaleExpansionZeroElim(B3, 2.0 * cdy,tmp1);
		scaleExpansionZeroElim(C1, cdytail,helper);
		scaleExpansionZeroElim(helper,bdx,tmp2);
		tmp2 = fastExpansionSumZeroElim(tmp1, tmp2);
		scaleExpansionZeroElim(C2, cdytail,helper);
		scaleExpansionZeroElim(helper, -adx,tmp1);
		tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
		F = fastExpansionSumZeroElim(F, tmp1);
	}
	if ((adxtail != 0.0) || (adytail != 0.0))
	{
		if ((bdxtail != 0.0) || (bdytail != 0.0) || (cdxtail != 0.0) || (cdytail != 0.0))
		{
			twoProduct(bdxtail, cdy, res1, err1);
			twoProduct(bdx, cdytail, res2, err2);
			A1[0] = res1;
			A1[1] = err1;
			A2[0] = res2;
			A2[1] = err2;
			tmp1 = twoTwoSum(A1, A2);
			twoProduct(cdxtail, -bdy, res1, err1);
			twoProduct(cdx, -bdytail, res2, err2);
			A1[0] = res1;
			A1[1] = err1;
			A2[0] = res2;
			A2[1] = err2;
			tmp2 = twoTwoSum(A1, A2);
			E1 = fastExpansionSumZeroElim(tmp1, tmp2);
			twoProduct(bdxtail, cdytail, res1, err1);
			twoProduct(cdxtail, bdytail, res2, err2);
			A1[0] = res1;
			A1[1] = err1;
			A2[0] = res2;
			A2[1] = err2;
			E2 = twoTwoDiff(A1, A2);
		}
		else
		{
			E1.push_back(0.0);
			E2.push_back(0.0);
		}
		if (adxtail != 0.0)
		{
			scaleExpansionZeroElim(D1, adxtail,tmp1);
			scaleExpansionZeroElim(E1, adxtail,E3);
			scaleExpansionZeroElim(E3, 2.0 * adx,tmp2);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
			if (bdytail != 0.0)
			{
				scaleExpansionZeroElim(C3, adxtail,tmp1);
				scaleExpansionZeroElim(tmp1, bdytail,tmp2);
				F = fastExpansionSumZeroElim(F, tmp2);
			}
			if (cdytail != 0.0)
			{
				scaleExpansionZeroElim(C2, -adxtail,tmp1);
				scaleExpansionZeroElim(tmp1, cdytail,tmp2);
				F = fastExpansionSumZeroElim(F, tmp2);
			}
			scaleExpansionZeroElim(E2, adxtail,D1);
			scaleExpansionZeroElim(D1, 2.0 * adx,tmp1);
			scaleExpansionZeroElim(D1, adxtail,tmp2);
			tmp2 = fastExpansionSumZeroElim(tmp1, tmp2);
			scaleExpansionZeroElim(E3, adxtail,tmp1);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
		}
		if (adytail != 0.0)
		{
			scaleExpansionZeroElim(D2, adytail,tmp1);
			scaleExpansionZeroElim(E1, adytail,D2);
			scaleExpansionZeroElim(D2, 2.0 * ady,tmp2);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
			scaleExpansionZeroElim(D2, adytail,tmp1);
			scaleExpansionZeroElim(E2, adytail,D1);
			scaleExpansionZeroElim(D1, 2.0 * ady,E1);
			scaleExpansionZeroElim(D1, adytail,E3);
			tmp2 = fastExpansionSumZeroElim(E1, E3);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
		}
	}
	if ((bdxtail != 0.0) || (bdytail != 0.0))
	{
		if ((cdxtail != 0.0) || (cdytail != 0.0) || (adxtail != 0.0) || (adytail != 0.0))
		{
			twoProduct(cdxtail, ady, res1, err1);
			twoProduct(cdx, adytail, res2, err2);
			A1[0] = res1;
			A1[1] = err1;
			A2[0] = res2;
			A2[1] = err2;
			tmp1 = twoTwoSum(A1, A2);
			twoProduct(adxtail, -cdy, res1, err1);
			twoProduct(adx, -cdytail, res2, err2);
			A1[0] = res1;
			A1[1] = err1;
			A2[0] = res2;
			A2[1] = err2;
			tmp2 = twoTwoSum(A1, A2);
			E1 = fastExpansionSumZeroElim(tmp1, tmp2);
			twoProduct(cdxtail, adytail, res1, err1);
			twoProduct(adxtail, cdytail, res2, err2);
			A1[0] = res1;
			A1[1] = err1;
			A2[0] = res2;
			A2[1] = err2;
			E2 = twoTwoDiff(A1, A2);
		}
		else
		{
			E1.push_back(0.0);
			E2.push_back(0.0);
		}
		if (bdxtail != 0.0)
		{
			scaleExpansionZeroElim(D3, bdxtail,tmp1);
			scaleExpansionZeroElim(E1, bdxtail,E3);
			scaleExpansionZeroElim(E3, 2.0 * bdx,tmp2);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
			if (cdytail != 0.0)
			{
				scaleExpansionZeroElim(C1, bdxtail,tmp1);
				scaleExpansionZeroElim(tmp1, cdytail,tmp2);
				F = fastExpansionSumZeroElim(F, tmp2);
			}
			if (adytail != 0.0)
			{
				scaleExpansionZeroElim(C3, -bdxtail,tmp1);
				scaleExpansionZeroElim(tmp1, adytail,tmp2);
				F = fastExpansionSumZeroElim(F, tmp2);
			}
			scaleExpansionZeroElim(E3, bdxtail,tmp1);
			scaleExpansionZeroElim(E2, bdxtail,D1);
			scaleExpansionZeroElim(D1, 2.0 * bdx,D2);
			scaleExpansionZeroElim(D1, bdxtail,D3);
			tmp2 = fastExpansionSumZeroElim(D2, D3);
			tmp2 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp2);
		}
		if (bdytail != 0.0)
		{
			scaleExpansionZeroElim(B1, bdytail,tmp1);
			scaleExpansionZeroElim(E1, bdytail,D1);
			scaleExpansionZeroElim(D1, 2.0 * bdy,tmp2);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
			scaleExpansionZeroElim(D1, bdytail,tmp1);
			scaleExpansionZeroElim(E2, bdytail,D1);
			scaleExpansionZeroElim(D1, 2.0 * bdy,D2);
			scaleExpansionZeroElim(D1, bdytail,D3);
			tmp2 = fastExpansionSumZeroElim(D2, D3);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
		}
	}
	if ((cdxtail != 0.0) || (cdytail != 0.0))
	{
		if ((adxtail != 0.0) || (adytail != 0.0) || (bdxtail != 0.0) || (bdytail != 0.0))
		{
			twoProduct(adxtail, bdy, res1, err1);
			twoProduct(adx, bdytail, res2, err2);
			A1[0] = res1;
			A1[1] = err1;
			A2[0] = res2;
			A2[1] = err2;
			tmp1 = twoTwoSum(A1, A2);
			twoProduct(bdxtail, -ady, res1, err1);
			twoProduct(bdx, -adytail, res2, err2);
			A1[0] = res1;
			A1[1] = err1;
			A2[0] = res2;
			A2[1] = err2;
			tmp2 = twoTwoSum(A1, A2);
			E1 = fastExpansionSumZeroElim(tmp1, tmp2);
			twoProduct(adxtail, bdytail, res1, err1);
			twoProduct(bdxtail, adytail, res2, err2);
			A1[0] = res1;
			A1[1] = err1;
			A2[0] = res2;
			A2[1] = err2;
			E2 = twoTwoDiff(A1, A2);
		}
		else
		{
			E1.push_back(0.0);
			E2.push_back(0.0);
		}
		if (cdxtail != 0.0)
		{
			scaleExpansionZeroElim(B2, cdxtail,tmp1);
			scaleExpansionZeroElim(E1, cdxtail,E3);
			scaleExpansionZeroElim(E3, 2.0 * cdx,tmp2);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
			if (adytail != 0.0)
			{
				scaleExpansionZeroElim(C2, cdxtail,tmp1);
				scaleExpansionZeroElim(tmp1, adytail,tmp2);
				F = fastExpansionSumZeroElim(F, tmp2);
			}
			if (bdytail != 0.0)
			{
				scaleExpansionZeroElim(C1, -cdxtail,tmp1);
				scaleExpansionZeroElim(tmp1, bdytail,tmp2);
				F = fastExpansionSumZeroElim(F, tmp2);
			}
			scaleExpansionZeroElim(E3, cdxtail,tmp1);
			scaleExpansionZeroElim(E2, cdxtail,D1);
			scaleExpansionZeroElim(D1, 2.0 * cdx,D2);
			scaleExpansionZeroElim(D1, cdxtail,D3);
			C1 = fastExpansionSumZeroElim(D2, D3);
			C2 = fastExpansionSumZeroElim(tmp1, C1);
			F = fastExpansionSumZeroElim(F, C2);
		}
		if (cdytail != 0.0)
		{
			scaleExpansionZeroElim(B3, cdytail,tmp1);
			scaleExpansionZeroElim(E1, cdytail,B1);
			scaleExpansionZeroElim(B1, 2.0 * cdy,tmp2);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
			scaleExpansionZeroElim(B1, cdytail,tmp1);
			scaleExpansionZeroElim(E2, cdytail,B1);
			scaleExpansionZeroElim(B1, 2.0 * cdy,B2);
			scaleExpansionZeroElim(B1, cdytail,B3);
			tmp2 = fastExpansionSumZeroElim(B2, B3);
			tmp1 = fastExpansionSumZeroElim(tmp1, tmp2);
			F = fastExpansionSumZeroElim(F, tmp1);
		}
	}
	return F.back();
}

double incircle(const Vector2D& point_1,
		const Vector2D& point_2,
		const Vector2D& point_3,
		const Vector2D& point_4)
{
	const double errBoundA = (10.0 + 96.0 * epsilon) * epsilon; // Acceptable error range for this calculation precision.

	const double adx = point_1.x - point_4.x;
	const double bdx = point_2.x - point_4.x;
	const double cdx = point_3.x - point_4.x;
	const double ady = point_1.y - point_4.y;
	const double bdy = point_2.y - point_4.y;
	const double cdy = point_3.y - point_4.y;

	const double alift = adx * adx + ady * ady;
	const double blift = bdx * bdx + bdy * bdy;
	const double clift = cdx * cdx + cdy * cdy;
	const double bdxcdy = bdx * cdy;
	const double cdxbdy = cdx * bdy;
	const double cdxady = cdx * ady;
	const double adxcdy = adx * cdy;
	const double adxbdy = adx * bdy;
	const double bdxady = bdx * ady;

	//Calculating regular 3x3 matrix determinant.
	const double det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);

	const double permanent = (fabs(bdxcdy) + fabs(cdxbdy)) * alift
	  + (fabs(cdxady) + fabs(adxcdy)) * blift
	  + (fabs(adxbdy) + fabs(bdxady)) * clift;
	const double errbound = errBoundA * permanent;
	if ((det > errbound) || (-det > errbound))
	{ //Determinant is out of error bounds (accurate enough).
		return det;
	}
	return incircleadapt(point_1,
			     point_2,
			     point_3,
			     point_4,
			     permanent); //calling the adaptive function.
}
