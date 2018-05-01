#ifndef DELAUNAY3D_HPP
#define DELAUNAY3D_HPP 1

#include "Vector3D.hpp"
#include "Tetrahedron.hpp"
#include <string>
#include <vector>
#include <set>
#include <array>
#include <boost/container/flat_set.hpp>

using std::vector;
using std::string;
//using std::stack;

class Delaunay3D
{
public:
	vector<Tetrahedron> tetras_;
	vector<Vector3D> points_;
	boost::container::flat_set<size_t> empty_tetras_;
	std::size_t Norg_;
	std::size_t outside_neighbor_;


	Delaunay3D();

	Delaunay3D(Delaunay3D const& other);

	Delaunay3D& operator=(Delaunay3D const& other);

	~Delaunay3D();

	void Build(vector<Vector3D> const& points,Vector3D const& maxv,Vector3D const& minv,std::vector<size_t> &order);

	void BuildExtra(vector<Vector3D> const& points);

	void output(string const& filename)const;

	bool CheckCorrect(void);

	void Clean(void);
private:
	void InsertPoint(std::size_t index);
	std::size_t Walk(std::size_t point, std::size_t first_guess);
	void flip14(std::size_t point,std::size_t tetra);
	void flip23(std::size_t tetra0, std::size_t tetra1,std::size_t location0,bool flat_check);
	void flip32(std::size_t tetra0, std::size_t tetra1, std::size_t location0, std::size_t shared_loction, bool flat_check);
	void flip44(std::size_t tetra0, std::size_t tetra1, std::size_t location0, std::size_t neigh0, std::size_t neigh1);
	void FindFlip(std::size_t tetra0, std::size_t tetra1, std::size_t p, size_t p_loc, size_t other_point_loc);
	void ExactFlip(std::size_t tetra0, std::size_t tetra1, std::size_t p);
	std::size_t FindThirdNeighbor(std::size_t tetra0, std::size_t tetra1);

	std::array<Vector3D, 3> b3_temp_,b3_temp2_;
	std::array<Vector3D, 4> b4_temp_;
	std::array<Vector3D, 5> b5_temp_;
	std::array<std::size_t, 4> b4s_temp_,b4s_temp2_;
	std::array<std::size_t, 8> b8s_temp_;
	vector<std::size_t> to_check_;
	std::size_t last_checked_;	
	Tetrahedron tet_temp0_, tet_temp1_, newtet_;
};

#endif //DELAUNAY3D_HPP