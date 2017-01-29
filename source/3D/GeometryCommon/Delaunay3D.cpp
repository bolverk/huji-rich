#include "Delaunay3D.hpp"
#include "Predicates3D.hpp"
#include <limits>
#include <algorithm>
#include <fstream>
#include "HilbertOrder3D.hpp"
#include <boost/foreach.hpp>

//#define runcheks 1

namespace
{
	Vector3D PlaneLineIntersection(boost::array<Vector3D, 3> const& plane, Vector3D const& A, Vector3D const& B)
	{
		Vector3D N = CrossProduct(plane[1] - plane[0], plane[2] - plane[0]);
		Vector3D mu = B - A;
		double m = ScalarProd(N, plane[0] - A) / ScalarProd(N, mu);
		return A + m*mu;
	}

	std::size_t GetOppositePoint(Tetrahedron const& tetra, std::size_t neighbor)
	{
		for(boost::array<size_t,4>::const_iterator it = tetra.neighbors.begin();it!= tetra.neighbors.end();++it)
			if (*it == neighbor)
				return static_cast<size_t>(it- tetra.neighbors.begin());
		assert(false);
	}

	std::size_t GetPointLocationInTetra(Tetrahedron const& tetra, std::size_t point)
	{
		for (std::size_t i = 0; i < 4; ++i)
			if (tetra.points[i] == point)
				return i;
		assert(false);
	}

	std::pair<std::size_t,double> InTriangle(boost::array<Vector3D, 3> const& triangle, Vector3D const& p)
	{
		// returns the smallest area of the cross product in units of the triangle area, negative if outside
		Vector3D a = CrossProduct(triangle[1] - triangle[0], p - triangle[0]);
		Vector3D b = CrossProduct(triangle[2] - triangle[1], p - triangle[1]);
		Vector3D c = CrossProduct(triangle[0] - triangle[2], p - triangle[2]);
		Vector3D d = CrossProduct(triangle[1] - triangle[0], triangle[2] - triangle[0]);
		double ad = ScalarProd(d, a);
		double bd = ScalarProd(d, b);
		double cd = ScalarProd(d, c);
		std::pair<std::size_t, double> res(0, 0);
		if (ad < 0)
			++res.first;
		if (bd < 0)
			++res.first;
		if (cd < 0)
			++res.first;
		res.second = std::min(std::min(std::abs(ad), std::abs(bd)), std::abs(cd)) / (ScalarProd(d, d));
		return res;
	}

	bool Are44(Tetrahedron const& T0, Tetrahedron const& T1, std::size_t loc_in_0, vector<Tetrahedron> const& tetras,
		std::size_t &N3,std::size_t &N4)
	{
		N3 = T0.neighbors[loc_in_0];
		N4 = T1.neighbors[GetPointLocationInTetra(T1, T0.points[loc_in_0])];
		for (std::size_t i = 0; i < 4; ++i)
			if (tetras[N3].neighbors[i] == N4)
				return true;
		return false;
	}

}
void Delaunay3D::flip23(std::size_t tetra0, std::size_t tetra1, std::size_t location0)
{
	bool used_empty = false;
	std::size_t Nloc = tetras_.size();
	if (!empty_tetras_.empty())
	{
		used_empty = true;
		Nloc = *empty_tetras_.begin();
		empty_tetras_.erase(empty_tetras_.begin());
	}

	// location is the location of the point in tetra that is not in the joint triangle
	Tetrahedron newtet, oldtet0(tetras_[tetra0]), oldtet1(tetras_[tetra1]);
	std::size_t location1 = GetOppositePoint(oldtet1, tetra0);

	newtet.points[0] = oldtet0.points[location0];
	newtet.points[1] = oldtet0.points[(location0 + 1)%4];
	newtet.points[2] = oldtet0.points[(location0 + 2) % 4];
	newtet.points[3] = oldtet1.points[location1];
	newtet.neighbors[0] = oldtet1.neighbors[GetPointLocationInTetra(
		oldtet1, oldtet0.points[(location0 + 3) % 4])];
	newtet.neighbors[1] = tetra0;
	newtet.neighbors[2] = tetra1;
	newtet.neighbors[3] = oldtet0.neighbors[(location0 + 3) % 4];
	if (newtet.neighbors[0] != outside_neighbor_)
		tetras_[newtet.neighbors[0]].neighbors[GetOppositePoint(tetras_[newtet.neighbors[0]], tetra1)] = Nloc;
	if (newtet.neighbors[3] != outside_neighbor_)
		tetras_[newtet.neighbors[3]].neighbors[GetOppositePoint(tetras_[newtet.neighbors[3]], tetra0)] = Nloc;
	if (used_empty)
		tetras_[Nloc] = newtet;
	else
		tetras_.push_back(newtet);

	tetras_[tetra0].points[0] =oldtet0.points[location0];
	tetras_[tetra0].points[3] = oldtet1.points[location1];
	tetras_[tetra0].points[1] = oldtet0.points[(location0 + 2) % 4];
	tetras_[tetra0].points[2] = oldtet0.points[(location0 + 3) % 4];
	tetras_[tetra0].neighbors[0] = oldtet1.neighbors[GetPointLocationInTetra(
		oldtet1, oldtet0.points[(location0 + 1) % 4])];
	tetras_[tetra0].neighbors[1] = tetra1;
	tetras_[tetra0].neighbors[2] = Nloc;
	tetras_[tetra0].neighbors[3] = oldtet0.neighbors[(location0 + 1) % 4];
	if (tetras_[tetra0].neighbors[0] != outside_neighbor_)
		tetras_[tetras_[tetra0].neighbors[0]].neighbors[GetOppositePoint(tetras_[tetras_[tetra0].neighbors[0]], tetra1)] = tetra0;

	tetras_[tetra1].points[0] = oldtet0.points[location0];
	tetras_[tetra1].points[3] = oldtet1.points[location1];
	tetras_[tetra1].points[1] = oldtet0.points[(location0 + 3) % 4];
	tetras_[tetra1].points[2] = oldtet0.points[(location0 + 1) % 4];
	tetras_[tetra1].neighbors[0] = oldtet1.neighbors[GetPointLocationInTetra(
		oldtet1, oldtet0.points[(location0 + 2) % 4])];
	tetras_[tetra1].neighbors[1] = Nloc;
	tetras_[tetra1].neighbors[2] = tetra0;
	tetras_[tetra1].neighbors[3] = oldtet0.neighbors[(location0 + 2) % 4];
	if (tetras_[tetra1].neighbors[3] != outside_neighbor_)
		tetras_[tetras_[tetra1].neighbors[3]].neighbors[GetOppositePoint(tetras_[tetras_[tetra1].neighbors[3]], tetra0)] = tetra1;

	for (int i = 0; i < 4; ++i)
		b4_temp_[static_cast<size_t>(i)] = points_[tetras_[Nloc].points[static_cast<size_t>(i)]];
	double orient = orient3d(b4_temp_);
	if (orient > 0)
	{
		std::size_t temp = tetras_[Nloc].points[0];
		tetras_[Nloc].points[0] = tetras_[Nloc].points[1];
		tetras_[Nloc].points[1] = temp;
		temp = tetras_[Nloc].neighbors[0];
		tetras_[Nloc].neighbors[0] = tetras_[Nloc].neighbors[1];
		tetras_[Nloc].neighbors[1] = temp;
	}
	else
	{
		if (std::abs(orient) == 0)
		{
			std::size_t my_loc = GetPointLocationInTetra(tetras_[Nloc], oldtet0.points[location0]);
			std::size_t neigh = tetras_[Nloc].neighbors[my_loc];
			std::size_t opp = GetOppositePoint(tetras_[neigh], Nloc);
			for (int i = 0; i < 3; ++i)
				b4_temp_[static_cast<size_t>(i)] = points_[tetras_[Nloc].points[(my_loc + static_cast<size_t>(i) + 1) % 4]];
			b4_temp_[3] = points_[tetras_[neigh].points[opp]];
			if (orient3d(b4_temp_) > 0)
			{
				std::size_t temp = tetras_[Nloc].points[0];
				tetras_[Nloc].points[0] = tetras_[Nloc].points[1];
				tetras_[Nloc].points[1] = temp;
				temp = tetras_[Nloc].neighbors[0];
				tetras_[Nloc].neighbors[0] = tetras_[Nloc].neighbors[1];
				tetras_[Nloc].neighbors[1] = temp;
			}
		}
	}


	for (int i = 0; i < 4; ++i)
		b4_temp_[i] = points_[tetras_[tetra0].points[static_cast<size_t>(i)]];
	orient = orient3d(b4_temp_);
	if (orient > 0)
	{
		std::size_t temp = tetras_[tetra0].points[0];
		tetras_[tetra0].points[0] = tetras_[tetra0].points[1];
		tetras_[tetra0].points[1] = temp;
		temp = tetras_[tetra0].neighbors[0];
		tetras_[tetra0].neighbors[0] = tetras_[tetra0].neighbors[1];
		tetras_[tetra0].neighbors[1] = temp;
	}
	else
	{
		if (std::abs(orient) == 0)
		{
			std::size_t my_loc = GetPointLocationInTetra(tetras_[tetra0], oldtet0.points[location0]);
			std::size_t neigh = tetras_[tetra0].neighbors[my_loc];
			std::size_t opp = GetOppositePoint(tetras_[neigh], tetra0);
			for (int i = 0; i < 3; ++i)
				b4_temp_[static_cast<size_t>(i)] = points_[tetras_[tetra0].points[(my_loc + static_cast<size_t>(i) + 1) % 4]];
			b4_temp_[3] = points_[tetras_[neigh].points[opp]];
			if (orient3d(b4_temp_) > 0)
			{
				std::size_t temp = tetras_[tetra0].points[0];
				tetras_[tetra0].points[0] = tetras_[tetra0].points[1];
				tetras_[tetra0].points[1] = temp;
				temp = tetras_[tetra0].neighbors[0];
				tetras_[tetra0].neighbors[0] = tetras_[tetra0].neighbors[1];
				tetras_[tetra0].neighbors[1] = temp;
			}
		}
	}


	for (int i = 0; i < 4; ++i)
		b4_temp_[static_cast<size_t>(i)] = points_[tetras_[tetra1].points[static_cast<size_t>(i)]];
	orient = orient3d(b4_temp_);
	if (orient > 0)
	{
		std::size_t temp = tetras_[tetra1].points[0];
		tetras_[tetra1].points[0] = tetras_[tetra1].points[1];
		tetras_[tetra1].points[1] = temp;
		temp = tetras_[tetra1].neighbors[0];
		tetras_[tetra1].neighbors[0] = tetras_[tetra1].neighbors[1];
		tetras_[tetra1].neighbors[1] = temp;
	}
	else
	{
		if (std::abs(orient) == 0)
		{
			std::size_t my_loc = GetPointLocationInTetra(tetras_[tetra1], oldtet0.points[location0]);
			std::size_t neigh = tetras_[tetra1].neighbors[my_loc];
			std::size_t opp = GetOppositePoint(tetras_[neigh], tetra1);
			for (int i = 0; i < 3; ++i)
				b4_temp_[static_cast<size_t>(i)] = points_[tetras_[tetra1].points[(my_loc + static_cast<size_t>(i) + 1) % 4]];
			b4_temp_[3] = points_[tetras_[neigh].points[opp]];
			if (orient3d(b4_temp_) > 0)
			{
				std::size_t temp = tetras_[tetra1].points[0];
				tetras_[tetra1].points[0] = tetras_[tetra1].points[1];
				tetras_[tetra1].points[1] = temp;
				temp = tetras_[tetra1].neighbors[0];
				tetras_[tetra1].neighbors[0] = tetras_[tetra1].neighbors[1];
				tetras_[tetra1].neighbors[1] = temp;
			}
		}
	}

	to_check_.push_back(tetra0);
	to_check_.push_back(tetra1);
	to_check_.push_back(Nloc);
}

void Delaunay3D::flip32(std::size_t tetra0, std::size_t tetra1, std::size_t location0,std::size_t shared_loction)
{
	// shared_loction is the point that is to be shared in the 2 new tetras. It is the point opposite to the shared edge by the three tetras in the triangle joint with the two tetras to check
	std::size_t other_point = 20, other_point2 = 20;
	for (std::size_t i = 0; i < 4; ++i)
	{
		if (i == location0 || i == shared_loction)
			continue;
		if (other_point == 20)
			other_point = i;
		else
			other_point2 = i;
	}

	std::size_t location1 = GetOppositePoint(tetras_[tetra1], tetra0);
	std::size_t third_tetra = tetras_[tetra0].neighbors[shared_loction];
	Tetrahedron newtet, old(tetras_[tetra0]);

	newtet.points[0] = tetras_[tetra0].points[location0];
	newtet.points[1] = tetras_[tetra0].points[shared_loction];
	newtet.points[2] = tetras_[tetra0].points[other_point];
	newtet.points[3] = tetras_[tetra1].points[location1];

	newtet.neighbors[0] = tetras_[tetra1].neighbors[GetPointLocationInTetra(
		tetras_[tetra1], tetras_[tetra0].points[other_point2])];
	newtet.neighbors[1] = tetras_[third_tetra].neighbors[GetPointLocationInTetra(
		tetras_[third_tetra], tetras_[tetra0].points[other_point2])];
	newtet.neighbors[2] = tetra1;
	newtet.neighbors[3] = tetras_[tetra0].neighbors[other_point2];
	if (newtet.neighbors[0] != outside_neighbor_)
		tetras_[newtet.neighbors[0]].neighbors[GetOppositePoint(tetras_[newtet.neighbors[0]], tetra1)] = tetra0;
	if (newtet.neighbors[1] != outside_neighbor_)
		tetras_[newtet.neighbors[1]].neighbors[GetOppositePoint(tetras_[newtet.neighbors[1]], third_tetra)] = tetra0;
	for (int i = 0; i < 4; ++i)
		b4_temp_[static_cast<size_t>(i)] = points_[newtet.points[static_cast<size_t>(i)]];
	if (orient3d(b4_temp_) > 0)
	{
		std::size_t temp = newtet.points[0];
		newtet.points[0] = newtet.points[1];
		newtet.points[1] = temp;
		temp = newtet.neighbors[0];
		newtet.neighbors[0] = newtet.neighbors[1];
		newtet.neighbors[1] = temp;
	}
	tetras_[tetra0] = newtet;
	
	newtet.points[0] = old.points[location0];
	newtet.points[1] = old.points[shared_loction];
	newtet.points[2] = old.points[other_point2];
	newtet.points[3] = tetras_[tetra1].points[location1];

	newtet.neighbors[0] = tetras_[tetra1].neighbors[GetPointLocationInTetra(
		tetras_[tetra1], old.points[other_point])];
	newtet.neighbors[1] = tetras_[third_tetra].neighbors[GetPointLocationInTetra(
		tetras_[third_tetra], old.points[other_point])];
	newtet.neighbors[2] = tetra0;
	newtet.neighbors[3] = old.neighbors[other_point];
	if (newtet.neighbors[3] != outside_neighbor_)
		tetras_[newtet.neighbors[3]].neighbors[GetOppositePoint(tetras_[newtet.neighbors[3]], tetra0)] = tetra1;
	if (newtet.neighbors[1] != outside_neighbor_)
		tetras_[newtet.neighbors[1]].neighbors[GetOppositePoint(tetras_[newtet.neighbors[1]], third_tetra)] = tetra1;
	for (int i = 0; i < 4; ++i)
		b4_temp_[static_cast<size_t>(i)] = points_[newtet.points[static_cast<size_t>(i)]];
	if (orient3d(b4_temp_) > 0)
	{
		std::size_t temp = newtet.points[0];
		newtet.points[0] = newtet.points[1];
		newtet.points[1] = temp;
		temp = newtet.neighbors[0];
		newtet.neighbors[0] = newtet.neighbors[1];
		newtet.neighbors[1] = temp;
	}
	tetras_[tetra1] = newtet;

	to_check_.push_back(tetra0);
	to_check_.push_back(tetra1);

	empty_tetras_.insert(third_tetra);
	if (third_tetra == last_checked_)
		last_checked_ = tetra0;
}

void Delaunay3D::flip44(std::size_t tetra0, std::size_t tetra1, std::size_t location0, std::size_t neigh0,std::size_t neigh1)
{
	std::size_t shared_location = 20;
	for (std::size_t i = 0; i < 4; ++i)
	{
		if (tetras_[neigh0].neighbors[i] == tetra0)
		{
			shared_location = i;
			break;
		}
	}
	assert(shared_location != 20);

	flip23(tetra0, tetra1, location0);
	location0 = 20;
	for (std::size_t i = 0; i < 4; ++i)
	{
		if (tetras_[neigh0].neighbors[i] == neigh1)
		{
			location0 = i;
			break;
		}
	}
	assert(location0 != 20);
	
	flip32(neigh0, neigh1, location0, shared_location);
}

Delaunay3D::Delaunay3D()
{
	to_check_.reserve(100);
}


Delaunay3D::~Delaunay3D()
{}

void Delaunay3D::BuildExtra(vector<Vector3D> const& points)
{
	vector<std::size_t> order = HilbertOrder3D(points);
	size_t Nstart = points_.size();
	points_.insert(points_.end(), points.begin(), points.end());

	assert(to_check_.empty());
	for (std::size_t i = 0; i < points.size(); ++i)
		InsertPoint(order[i] + Nstart);
}

void Delaunay3D::Build(vector<Vector3D> const & points, Vector3D const& maxv, Vector3D const& minv)
{
	std::size_t Norg = points.size();
	Norg_ = Norg;
	points_.reserve(static_cast<std::size_t>(Norg+std::pow(Norg,0.6666)*14));
//	points_ = points;
	points_.assign(points.begin(), points.end());
	// Create large tetra points
	double factor = 50;
	points_.push_back(Vector3D(minv.x - 1.01*factor * (maxv.x - minv.x), minv.y - factor * (maxv.y - minv.y), minv.z - factor * (maxv.z - minv.z)));
	points_.push_back(Vector3D(0.5*(minv.x + maxv.x), maxv.y + 1.02*factor *(maxv.y - minv.y), minv.z - factor * (maxv.z - minv.z)));
	points_.push_back(Vector3D(maxv.x + 0.99*factor * (maxv.x - minv.x), minv.y - factor * (maxv.y - minv.y), minv.z - factor * (maxv.z - minv.z)));
	points_.push_back(Vector3D(0.5*(minv.x + maxv.x), 0.5*(minv.y + maxv.y), maxv.z + factor * (maxv.z - minv.z)));
	// Create large tetra
	outside_neighbor_=std::numeric_limits<std::size_t>::max();
	Tetrahedron tetra;
	tetra.points[0] = Norg;
	tetra.points[1] = Norg+2;
	tetra.points[2] = Norg+1;
	tetra.points[3] = Norg+3;
	tetra.neighbors[0] = outside_neighbor_;
	tetra.neighbors[1] = outside_neighbor_;
	tetra.neighbors[2] = outside_neighbor_;
	tetra.neighbors[3] = outside_neighbor_;
	tetras_.reserve(points_.capacity() * 7);
	tetras_.push_back(tetra);
	last_checked_ = 0;
	vector<std::size_t> order = HilbertOrder3D(points);
	
	assert(to_check_.empty());
	for (std::size_t i = 0; i < Norg; ++i)
		InsertPoint(order[i]);
}

void Delaunay3D::output(string const & filename) const
{
	std::ofstream fh(filename.c_str(), std::ostream::binary);

	std::size_t temp = tetras_.size() - empty_tetras_.size();
	fh.write(reinterpret_cast<const char*>(&temp), sizeof(std::size_t));
	temp = points_.size();
	fh.write(reinterpret_cast<const char*>(&Norg_), sizeof(std::size_t));
	fh.write(reinterpret_cast<const char*>(&temp), sizeof(std::size_t));

	for (std::size_t i = 0; i<points_.size(); ++i) 
	{
		fh.write(reinterpret_cast<const char*>(&points_[i].x), sizeof(double));
		fh.write(reinterpret_cast<const char*>(&points_[i].y), sizeof(double));
		fh.write(reinterpret_cast<const char*>(&points_[i].z), sizeof(double));
	}

	for (std::size_t i = 0; i<tetras_.size(); ++i) 
	{
		if (empty_tetras_.find(i) == empty_tetras_.end())
		{
			for (std::size_t j = 0; j < 4; ++j)
			{
				fh.write(reinterpret_cast<const char*>(&tetras_[i].points[j]), sizeof(std::size_t));
			}
		}
	}

	fh.close();
}

std::size_t Delaunay3D::FindThirdNeighbor(std::size_t tetra0,std::size_t tetra1)
{
	for (size_t i = 0; i < 4; ++i)
	{
		for (size_t j = 0; j < 4; ++j)
		{
			if (tetras_[tetra0].neighbors[i] == tetras_[tetra1].neighbors[j])
			{
				b8s_temp_[0] = tetras_[tetra1].neighbors[j];
				return 1;
			}
		}
	}
	return 0;
}

void Delaunay3D::ExactFlip(std::size_t tetra0, std::size_t tetra1, std::size_t p)
{
	std::size_t out_counter = 0;
	std::size_t other_point = 0;
	std::size_t in_counter = 0;
	std::size_t flat_counter = 0;
	boost::array<std::size_t, 3> out_check;
	Tetrahedron const& T0 = tetras_[tetra0];
	Tetrahedron const& T1 = tetras_[tetra1];
	std::size_t p_loc = GetPointLocationInTetra(T0, p);
	for (int i = 0; i < 3; ++i)
		b3_temp_[i] = points_[T0.points[(p_loc + 1 + static_cast<size_t>(i)) % 4]];

	for (std::size_t i = 0; i < 4; ++i)
	{
		if (T1.neighbors[i] == tetra0)
		{
			other_point = T1.points[i];
			break;
		}
	}

	b4_temp_[0] = b3_temp_[1];
	b4_temp_[1] = b3_temp_[2];
	b4_temp_[2] = points_[p];
	b4_temp_[3] = b3_temp_[0];
	double test0 = orient3d(b4_temp_);
	if (std::abs(test0) == 0)
		in_counter++;
	b4_temp_[3] = points_[other_point];
	double test1 = orient3d(b4_temp_);
	if (std::abs(test1) == 0)
		flat_counter++;
	if (test0*test1 > 0)
		out_check[0] = 0;
	else
	{
		out_check[0] = 1;
		++out_counter;
	}	

	b4_temp_[0] = b3_temp_[2];
	b4_temp_[1] = b3_temp_[0];
	b4_temp_[2] = points_[p];
	b4_temp_[3] = b3_temp_[1];
	test0 = orient3d(b4_temp_);
	if (std::abs(test0) == 0)
		in_counter++;
	b4_temp_[3] = points_[other_point];
	test1 = orient3d(b4_temp_);
	if (std::abs(test1) == 0)
		flat_counter++;
	if (test0*test1 > 0)
		out_check[1] = 0;
	else
	{
		out_check[1] = 1;
		++out_counter;
	}

	b4_temp_[0] = b3_temp_[0];
	b4_temp_[1] = b3_temp_[1];
	b4_temp_[2] = points_[p];
	b4_temp_[3] = b3_temp_[2];
	test0 = orient3d(b4_temp_);
	if (std::abs(test0) == 0)
		in_counter++;
	b4_temp_[3] = points_[other_point];
	test1 = orient3d(b4_temp_);
	if (std::abs(test1) == 0)
		flat_counter++;
	if (test0*test1 > 0)
		out_check[2] = 0;
	else
	{
		out_check[2] = 1;
		++out_counter;
	}
		
	if (out_counter == 0)
		flip23(tetra0, tetra1, p_loc);
	else
	{
		if (out_counter == 1)
		{
			std::size_t N_shared = FindThirdNeighbor(tetra0, tetra1);
			if (N_shared == 1)
			{
				std::size_t shared_loc = GetOppositePoint(tetras_[tetra0], b8s_temp_[0]);
				flip32(tetra0, tetra1, p_loc, shared_loc);
				return;
			}
			else
			{
				if (flat_counter > 0)
				{
					for (std::size_t i = 0; i < 3; ++i)
					{
						std::size_t N3, N4;
						if (out_check[i] == 1 && Are44(T0, T1, (p_loc + i + 1) % 4, tetras_, N3, N4))
						{
							flip44(tetra0, tetra1, p_loc, N3, N4);
							return;
						}
					}
				}
			}
		}
		else
		{
			if (in_counter == 3) // tetra0 is flat
			{
				//is a 32 flip possible?
				std::size_t N_shared = FindThirdNeighbor(tetra0, tetra1);
				if (N_shared == 1)
				{
					std::size_t shared_loc = GetOppositePoint(tetras_[tetra0], b8s_temp_[0]);
					flip32(tetra0, tetra1, p_loc, shared_loc);
					return;
				}
				else
					flip23(tetra0, tetra1, p_loc);
			}
		}
	}
}

void Delaunay3D::FindFlip(std::size_t tetra0,std::size_t tetra1,std::size_t p)
{
	std::size_t p_loc = GetPointLocationInTetra(tetras_[tetra0], p);
	std::size_t other_point_loc = GetOppositePoint(tetras_[tetra1], tetra0);
	for (std::size_t i = 0; i < 3;++i)
		b3_temp_[i] = points_[tetras_[tetra0].points[(p_loc + i + 1) % 4]];
	Vector3D intersection = PlaneLineIntersection(b3_temp_, points_[p], 
		points_[tetras_[tetra1].points[other_point_loc]]);
	std::pair<std::size_t, double> outside_intersection = InTriangle(b3_temp_, intersection);
	
	if (outside_intersection.second<1e-6)
	{
		ExactFlip(tetra0, tetra1, p);
		return;
	}
	if (outside_intersection.first == 0)
	{
		flip23(tetra0, tetra1, p_loc);
		return;
	}
	if (outside_intersection.first == 1)
	{
		// Do we have a shared neighbor?
		std::size_t Nshared = FindThirdNeighbor(tetra0, tetra1);
		if (Nshared == 1)
		{
			std::size_t shared_loc = GetOppositePoint(tetras_[tetra0], b8s_temp_[0]);
			flip32(tetra0, tetra1, p_loc, shared_loc);
		}
	}	
}

void Delaunay3D::InsertPoint(std::size_t index)
{
	std::size_t to_split = Walk(index, last_checked_);
	last_checked_ = to_split;
	flip14(index, to_split);
	while (!to_check_.empty())
	{
		std::size_t cur_check = to_check_.back();
		to_check_.pop_back();
		std::size_t to_flip = tetras_[cur_check].neighbors[GetPointLocationInTetra(tetras_[cur_check], index)];
		if (to_flip == outside_neighbor_ || (empty_tetras_.find(cur_check) != empty_tetras_.end()))
			continue;
		/// check here that we have correct sign for insphere test !!
		b5_temp_[0] = points_[tetras_[cur_check].points[0]];
		b5_temp_[1] = points_[tetras_[cur_check].points[1]];
		b5_temp_[2] = points_[tetras_[cur_check].points[2]];
		b5_temp_[3] = points_[tetras_[cur_check].points[3]];
		b5_temp_[4] = points_[tetras_[to_flip].points[GetOppositePoint(tetras_[to_flip], cur_check)]];
		if (insphere(b5_temp_) < 0)
		{
			FindFlip(cur_check, to_flip,index);
		}
	}
}

std::size_t Delaunay3D::Walk(std::size_t point, std::size_t first_guess) 
{
	bool good = false;
	std::size_t cur_facet = first_guess;
	std::size_t counter=0;
	b4_temp_[3] = points_[point];
	while (!good)
	{
		++counter;
		good = true;
		for (std::size_t i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 3; ++j)
				b4_temp_[static_cast<size_t>(j)] = points_[tetras_[cur_facet].points[(i + static_cast<size_t>(j) + 1) % 4]];
			int sign = 2 * static_cast<int>(i % 2) - 1;
			if ((orient3d(b4_temp_)*sign)>0)
			{
				good = false;
				cur_facet = tetras_[cur_facet].neighbors[i];
				break;
			}
		}
		assert(counter < 1e7);
	}
	return cur_facet;
}

void Delaunay3D::flip14(std::size_t point, std::size_t tetra)
{
	Tetrahedron toadd;
	boost::array<std::size_t, 3> Nloc;
	bool cleared_empty = false;
	if (empty_tetras_.size()>3)
	{
		cleared_empty = true;
		for (int i = 0; i < 3; ++i)
		{
			Nloc[static_cast<size_t>(i)] = *empty_tetras_.begin();
			empty_tetras_.erase(empty_tetras_.begin());
		}	
	}
	else
	{
		size_t Ntet = tetras_.size();
		for (int i = 0; i < 3; ++i)
			Nloc[static_cast<size_t>(i)] = Ntet + i;
	}

	toadd.neighbors[0] = Nloc[1];
	toadd.neighbors[1] = Nloc[2];
	toadd.neighbors[2] = tetras_[tetra].neighbors[0];
	toadd.neighbors[3] = tetra;
	toadd.points[0] = tetras_[tetra].points[1];
	toadd.points[1] = tetras_[tetra].points[2];
	toadd.points[2] = point;
	toadd.points[3] = tetras_[tetra].points[3];
	if(toadd.neighbors[2]!=outside_neighbor_)
		tetras_[toadd.neighbors[2]].neighbors[GetOppositePoint(tetras_[toadd.neighbors[2]], tetra)] = Nloc[0];
	if (!cleared_empty)
		tetras_.push_back(toadd);
	else
		tetras_[Nloc[0]] = toadd;

	toadd.neighbors[0] = Nloc[0];
	toadd.neighbors[1] = Nloc[2];
	toadd.neighbors[2] = tetra;
	toadd.neighbors[3] = tetras_[tetra].neighbors[1];
	toadd.points[0] = tetras_[tetra].points[0];
	toadd.points[1] = tetras_[tetra].points[2];
	toadd.points[2] = tetras_[tetra].points[3];
	toadd.points[3] = point;
	if (toadd.neighbors[3] != outside_neighbor_)
		tetras_[toadd.neighbors[3]].neighbors[GetOppositePoint(tetras_[toadd.neighbors[3]], tetra)] = Nloc[1];
	if (!cleared_empty)
		tetras_.push_back(toadd);
	else
		tetras_[Nloc[1]] = toadd;

	toadd.neighbors[0] = Nloc[0];
	toadd.neighbors[1] = tetra;
	toadd.neighbors[2] = Nloc[1];
	toadd.neighbors[3] = tetras_[tetra].neighbors[2];
	toadd.points[0] = tetras_[tetra].points[0];
	toadd.points[1] = tetras_[tetra].points[3];
	toadd.points[2] = tetras_[tetra].points[1];
	toadd.points[3] = point;
	if (toadd.neighbors[3] != outside_neighbor_)
		tetras_[toadd.neighbors[3]].neighbors[GetOppositePoint(tetras_[toadd.neighbors[3]], tetra)] = Nloc[2];
	if (!cleared_empty)
		tetras_.push_back(toadd);
	else
		tetras_[Nloc[2]] = toadd;


	tetras_[tetra].points[3] = point;
	tetras_[tetra].neighbors[0] = Nloc[0];
	tetras_[tetra].neighbors[1] = Nloc[1];
	tetras_[tetra].neighbors[2] = Nloc[2];
	
	to_check_.push_back(tetra);
	to_check_.push_back(Nloc[0]);
	to_check_.push_back(Nloc[1]);
	to_check_.push_back(Nloc[2]);
#ifdef runcheks
	for (std::size_t i = 0; i < 4; ++i)
		b4_temp_[i] = points_[tetras_[tetras_.size() - 1].points[i]];
	assert(orient3d(b4_temp_) <= 0);
	for (std::size_t i = 0; i < 4; ++i)
		b4_temp_[i] = points_[tetras_[tetras_.size()-2].points[i]];
	assert(orient3d(b4_temp_) <= 0);
	for (std::size_t i = 0; i < 4; ++i)
		b4_temp_[i] = points_[tetras_[tetras_.size() - 3].points[i]];
	assert(orient3d(b4_temp_) <= 0);
	for (std::size_t i = 0; i < 4; ++i)
		b4_temp_[i] = points_[tetras_[tetra].points[i]];
	assert(orient3d(b4_temp_) <= 0);
#endif
}

bool Delaunay3D::CheckCorrect(void)
{
	std::size_t Ntetra = tetras_.size();
	
	for (std::size_t i = 0; i < Ntetra; ++i)
	{
		if (empty_tetras_.find(i) != empty_tetras_.end())
			continue;
		Tetrahedron const& T = tetras_[i];		
		b5_temp_[0] = points_[T.points[0]];
		b5_temp_[1] = points_[T.points[1]];
		b5_temp_[2] = points_[T.points[2]];
		b5_temp_[3] = points_[T.points[3]];
		for (std::size_t j = 0; j < 4; ++j)
		{
			// Check same neighbors
			if (T.neighbors[j] != outside_neighbor_)
				GetOppositePoint(tetras_[T.neighbors[j]], i);
			// Check insphere
			if (T.neighbors[j] != outside_neighbor_)
			{
				std::size_t other = GetOppositePoint(tetras_[T.neighbors[j]], i);
				b5_temp_[4] = points_[tetras_[T.neighbors[j]].points[other]];
				assert(!(insphere(b5_temp_) < 0));
			}
		}			
	}
	return true;
}

void Delaunay3D::Clean(void)
{
	tetras_.clear();
	points_.clear();
	empty_tetras_.clear();
}
