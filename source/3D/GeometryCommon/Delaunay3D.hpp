#ifndef DELAUNAY3D_HPP
#define DELAUNAY3D_HPP 1

#include "Vector3D.hpp"
#include "Tetrahedron.hpp"
#include <string>
#include <vector>
#include <stack>
#include <set>


using std::vector;
using std::string;
using std::stack;

class Delaunay3D
{
public:
	vector<Tetrahedron> tetras_;
	vector<Vector3D> points_;
	std::set<std::size_t> empty_tetras_;
	std::size_t Norg_;
	std::size_t outside_neighbor_;


	Delaunay3D();

	~Delaunay3D();

	void Build(vector<Vector3D> const& points,Vector3D const& maxv,Vector3D const& minv);

	void BuildExtra(vector<Vector3D> const& points);

	void output(string const& filename)const;

	bool CheckCorrect(void);

	void Clean(void);
private:
	void InsertPoint(std::size_t index);
	std::size_t Walk(std::size_t point, std::size_t first_guess);
	void flip14(std::size_t point,std::size_t tetra);
	void flip23(std::size_t tetra0, std::size_t tetra1,std::size_t location0);
	void flip32(std::size_t tetra0, std::size_t tetra1, std::size_t location0, std::size_t shared_loction);
	void flip44(std::size_t tetra0, std::size_t tetra1, std::size_t location0, std::size_t neigh0, std::size_t neigh1);
	void FindFlip(std::size_t tetrao, std::size_t tetra1, std::size_t p);
	void ExactFlip(std::size_t tetra0, std::size_t tetra1, std::size_t p);
	std::size_t FindThirdNeighbor(std::size_t tetra0, std::size_t tetra1);

	boost::array<Vector3D, 3> b3_temp_,b3_temp2_;
	boost::array<Vector3D, 4> b4_temp_;
	boost::array<Vector3D, 5> b5_temp_;
	boost::array<std::size_t, 4> b4s_temp_,b4s_temp2_;
	boost::array<std::size_t, 8> b8s_temp_;
	stack<std::size_t> to_check_;
	std::size_t last_checked_;
};

#endif //DELAUNAY3D_HPP