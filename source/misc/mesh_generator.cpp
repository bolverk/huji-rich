#include "mesh_generator.hpp"

using namespace std;
typedef boost::mt19937_64 base_generator_type;

std::vector<Vector2D> RandSquare(int PointNum, Vector2D const& lowerleft, Vector2D const& upperright)
{
	return RandSquare(PointNum, lowerleft.x, upperright.x, lowerleft.y, upperright.y);
}

vector<Vector2D> RandPointsRmax(int PointNum, double Rmin, double Rmax,
	double xc, double yc)
{
	double ran[2];
	int counter = 0;
	vector<Vector2D> res;
	Vector2D point;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	while (counter < PointNum)
	{
		ran[0] = dist(generator);
		ran[1] = dist(generator);
		point.x = ((Rmax - Rmin)*ran[0] + Rmin)*cos(ran[1] * 2 * M_PI) + xc;
		point.y = ((Rmax - Rmin)*ran[0] + Rmin)*sin(ran[1] * 2 * M_PI) + yc;
		res.push_back(point);
		++counter;
	}
	return res;
}

vector<Vector2D> RandSquare(int PointNum, double xl, double xr, double yd,
	double yu)
{
	double ran[2];
	double xc = xr - xl;
	double yc = yu - yd;
	vector<Vector2D> res;
	Vector2D point;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	for (int i = 0; i < PointNum; ++i)
	{
		ran[0] = dist(generator);
		ran[1] = dist(generator);
		point.x = ran[0] * xc + xl;
		point.y = ran[1] * yc + yd;
		res.push_back(point);
	}
	return res;
}

std::vector<Vector2D> RandPointsCylinder(int PointNum, Vector2D const& ll, Vector2D const& ur)
{
	base_generator_type generator;
	double ran[2];
	double dx = ur.x - ll.x;
	double A = 1.0 / log(ur.y / ll.y);
	vector<Vector2D> res(static_cast<size_t>(PointNum));
	boost::random::uniform_real_distribution<> dist;
	for (int i = 0; i < PointNum; ++i)
	{
		ran[0] = dist(generator);
		ran[1] = dist(generator);
		res[i] = Vector2D(ran[0] * dx - ll.x, exp(ran[1]/A)*ll.y);
	}
	return res;
}

vector<Vector2D> RandSquare(int PointNum, boost::random::mt19937 &eng, double xl, double xr, double yd,
	double yu)
{
	double ran[2];
	double xc = xr - xl;
	double yc = yu - yd;
	vector<Vector2D> res;
	Vector2D point;
	boost::random::uniform_real_distribution<> dist;
	for (int i = 0; i < PointNum; ++i)
	{
		ran[0] = dist(eng);
		ran[1] = dist(eng);
		point.x = ran[0] * xc + xl;
		point.y = ran[1] * yc + yd;
		res.push_back(point);
	}
	return res;
}

vector<Vector2D> RandPointsRa(int PointNum, double Rmin, double Rmax, double alpha,
	Vector2D const& lowerleft, Vector2D const& upperright, Vector2D const& center)
{
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	double ran[2];
	double A = pow(Rmax*Rmin, alpha)*(alpha - 1) / (pow(Rmax, alpha)*Rmin - Rmax*pow(Rmin, alpha));
	int counter = 0;
	vector<Vector2D> res;
	Vector2D point;
	res.reserve(static_cast<size_t>(PointNum));
	while (counter < PointNum)
	{
		ran[0] = pow(pow(Rmin, 1 - alpha) - dist(generator)*(alpha - 1) / A, 1.0 / (1 - alpha));
		ran[1] = dist(generator);
		point.x = ran[0] * cos(ran[1] * 2 * M_PI) + center.x;
		point.y = ran[0] * sin(ran[1] * 2 * M_PI) + center.y;
		if ((point.x<upperright.x) && (point.x>lowerleft.x) && (point.y < upperright.y)
			&& (point.y > lowerleft.y))
		{
			res.push_back(point);
			++counter;
		}
	}
	return res;
}

vector<Vector2D> RandPointsR(int PointNum, double xl, double xr, double yd,
	double yu, double minR, double xc, double yc)
{
	double ran[2];
	int counter = 0;
	double maxR = sqrt(2.0)*max(max(max(xr - xc, yu - yc), yc - yd), xc - xl);
	vector<Vector2D> res;
	Vector2D point;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	while (counter < PointNum)
	{
		ran[0] = dist(generator);
		ran[1] = dist(generator);
		point.x = maxR*ran[0] * cos(ran[1] * 2 * M_PI) + xc;
		point.y = maxR*ran[0] * sin(ran[1] * 2 * M_PI) + yc;
		if (abs(point) > minR)
		{
			if (point.x > xl&&point.x<xr)
				if (point.y>yd&&point.y < yu)
				{
					res.push_back(point);
					++counter;
				}
		}
	}
	return res;
}

vector<Vector2D> CirclePointsRmax_1(int PointNum, double Rmin, double Rmax,
	double xc, double yc, double xmax, double ymax, double xmin, double ymin)
{
	const double Nphid = sqrt(PointNum * 2 * M_PI / log(Rmax / Rmin));
	const int Nr = static_cast<int>(Nphid*log(Rmax / Rmin) / (2 * M_PI));
	const int Nphi = static_cast<int>(floor(Nphid + 0.5));
	const double dphi = 2 * M_PI / Nphi;
	vector<Vector2D> res;
	for (int i = 0; i < Nr; ++i)
	{
		const double r = Rmin*exp(2 * M_PI*i / Nphi);
		for (int j = 0; j < Nphi; ++j)
		{
			const Vector2D pos(r*cos(dphi*j) + xc, r*sin(dphi*j) + yc);
			if (pos.x<xmax&&pos.x>xmin && pos.y<ymax&&pos.y>ymin)
				res.push_back(pos);
		}
	}
	return res;
}

vector<Vector2D> CirclePointsRmax_2(int PointNum, double Rmin, double Rmax,
	double xc, double yc, double xmax, double ymax, double xmin, double ymin)
{
	const double x = sqrt(Rmin / Rmax);
	const int Nr = int(sqrt(2 * PointNum*(1 - x) / (1 + x) / M_PI));
	const double A = 4 * pow(1 / sqrt(Rmin) - 1 / sqrt(Rmax), 2) / (Nr*Nr);
	vector<Vector2D> res;
	for (int i = 0; i < Nr; ++i)
	{
		const double r = 4 * Rmin / pow(2 - sqrt(A*Rmin)*i, 2);
		const double Nphi = int(2 * M_PI / sqrt(A*r));
		const double dphi = double(2 * M_PI / Nphi);
		for (int j = 0; j < Nphi; ++j)
		{
			const Vector2D pos(r*cos(dphi*j) + xc, r*sin(dphi*j) + yc);
			if (pos.x<xmax&&pos.x>xmin && pos.y<ymax&&pos.y>ymin)
				res.push_back(pos);
		}
	}
	return res;
}

vector<Vector2D> CirclePointsRmax(int PointNum, double Rmin, double Rmax,
	double xmax, double ymax, double xmin, double ymin, double xc, double yc)
{
	const double A = sqrt(M_PI*(Rmax*Rmax - Rmin*Rmin) / PointNum);
	const int Nr = int((Rmax - Rmin) / A);
	const double dr = (Rmax - Rmin) / Nr;
	vector<Vector2D> res;
	for (int i = 0; i < Nr; ++i)
	{
		const double r = Rmin + i*dr;
		const int Nphi = int(2 * M_PI*r / A);
		const double dphi = A / r;
		for (int j = 0; j < Nphi; ++j)
		{
			const Vector2D pos(r*cos(dphi*j) + xc, r*sin(dphi*j) + yc);
			if (pos.x<xmax&&pos.x>xmin && pos.y<ymax&&pos.y>ymin)
				res.push_back(pos);
		}
	}
	return res;
}

vector<Vector2D> cartesian_mesh(int nx, int ny,
	Vector2D const& lower_left,
	Vector2D const& upper_right)
{
	assert(nx > 0);
	assert(ny > 0);
	assert(upper_right.x > lower_left.x);
	assert(upper_right.y > lower_left.y);

	vector<Vector2D> res;
	const double dx = (upper_right.x - lower_left.x) /
		static_cast<double>(nx);
	const double dy = (upper_right.y - lower_left.y) /
		static_cast<double>(ny);
	for (double x = lower_left.x + 0.5*dx; x < upper_right.x; x += dx){
		for (double y = lower_left.y + 0.5*dy; y < upper_right.y; y += dy)
			res.push_back(Vector2D(x, y));
	}
	return res;
}

vector<Vector2D> circle_circumference
(size_t point_number,
double radius,
Vector2D const& center)
{
	vector<Vector2D> res(point_number);
	for (size_t i = 0; i < point_number; ++i)
	{
		const double angle = 2 * M_PI*double(i) / double(point_number);
		res[i] = center + pol2cart(radius, angle);
	}
	return res;
}

vector<Vector2D> Line(int PointNum, double xmin, double xmax, double ymin, double ymax)
{
	double dy = ymax - ymin, dx = xmax - xmin;
	double angle = atan2(dy, dx);
	double length = sqrt(dy*dy + dx*dx);
	double dl = length / PointNum;

	vector<Vector2D> res;
	Vector2D temp;

	for (int i = 0; i < PointNum; ++i)
	{
		temp.Set(xmin + dl*(i + 0.5)*cos(angle), ymin + dl*(i + 0.5)*sin(angle));
		res.push_back(temp);
	}
	return res;
}

vector<Vector2D> CirclePointsRmax_a(int PointNum,
	double Rmin,
	double Rmax,
	double xc,
	double yc,
	double xmax,
	double ymax,
	double xmin,
	double ymin,
	double alpha)
{
	const double N0 = sqrt(PointNum * 4 * M_PI*(alpha + 1) /
		(pow(Rmax, 2 * (alpha + 1)) -
		pow(Rmin, 2 * (alpha + 1))));
	const int Nr = int(floor((pow(Rmax, alpha + 1) -
		pow(Rmin, alpha + 1))*N0 / (2 * M_PI*(alpha + 1))));
	vector<Vector2D> res;
	for (int i = 0; i < Nr; ++i)
	{
		const double r = pow(2 * M_PI*i*(alpha + 1) / N0 +
			pow(Rmin, alpha + 1), 1.0 / (alpha + 1));
		const int Nphi = int(floor(N0*pow(r, 1 + alpha) + 1.5));
		const double dphi = 2 * M_PI / Nphi;
		for (int j = 0; j < Nphi; ++j)
		{
			const Vector2D pos(r*cos(dphi*j) + xc, r*sin(dphi*j) + yc);
			if (pos.x<xmax&&pos.x>xmin && pos.y<ymax&&pos.y>ymin)
				res.push_back(pos);
		}
	}
	return res;
}
