#include "find_affected_cells.hpp"
#include "geometry.hpp"
#include "boost/foreach.hpp"
#include "../misc/utils.hpp"
#include <cmath>

namespace
{
	bool periodic_duplicate_before(std::vector<Vector2D> const& toadd, Vector2D const& tocheck,Vector2D const& ll,Vector2D const& ur)
	{
		size_t N = toadd.size();
		if (N == 0)
			return false;
		if (tocheck.x > ur.x || tocheck.x<ll.x || tocheck.y>ur.y || tocheck.y < ll.y)
			return true;
		double dx = 1e-20;
		for (size_t i = 0; i < N; ++i)
			dx = std::max(dx, fastabs(toadd[i]));
		for (size_t i = 0; i < N; ++i)
		{
			double dV = fastabs(toadd[i] - tocheck);
			if (dV < dx*1e-6)
			{
				return true;
			}
		}
		return false;
	}

	bool contains(const vector<int>& v,	const int n)
	{
		for (size_t i = 0; i < v.size(); ++i) 
		{
			if (n == v[i])
				return true;
		}
		return false;
	}

	bool contains(const vector<int>& v, const std::vector<Vector2D> & added,Vector2D const& toadd, const int n,
		Vector2D const&ll, Vector2D const& ur)
	{
		if (fabs(toadd.x) > (1.01*(ur.x - ll.x)) || fabs(toadd.y) > (1.01*(ur.y - ll.y)))
			return true;
		for (size_t i = 0; i < v.size(); ++i)
		{
			if (n == v[i])
			{
				if(fastabs(added[i]-toadd)<1e-6*std::max(1e-20,std::max(fastabs(toadd),fastabs(added[i]))))
					return true;
			}
		}
		return false;
	}

	/*double bracket(double low, double num, double high)
	{
		return std::max(std::min(num, high), low);
	}*/
}

bool edge_circle_intersect
(const Edge& edge,
	const Circle& circle)
{
	Vector2D norm = normalize(edge.vertices.second - edge.vertices.first);
	Vector2D tocenter = circle.getCenter() - edge.vertices.first;
	Vector2D perpindicular = tocenter - ScalarProd(norm, tocenter) * norm;
	return ScalarProd(perpindicular, perpindicular) < circle.getRadius() * circle.getRadius();
}

namespace {
	bool cell_circle_intersect(const Tessellation& tess,
		int index,
		const Circle& circle)
	{
		BOOST_FOREACH(int ei, tess.GetCellEdges(index))
		{
			if (edge_circle_intersect(tess.GetEdge(ei), circle))
				return true;
		}
		return false;
	}

	void find_affected_cells2(const Tessellation& tess, int index, Circle circle,
		vector<int> &res, vector<int>& visited, std::vector<Vector2D> &added, bool periodic, Vector2D const& toaddsingle,
		std::vector<Vector2D>& visited_add,Vector2D const& ll,Vector2D const& ur)
	{
		if (periodic)
		{
			if (contains(visited, visited_add, toaddsingle, index,ll,ur))
				return;
		}
		else
			if (contains(visited, index))
				return;
		visited.push_back(index);
		visited_add.push_back(toaddsingle);
		if (!cell_circle_intersect(tess, tess.GetOriginalIndex(index), circle))
			return;
		else
		{
			res.push_back(index);
			if (periodic)
					added.push_back(toaddsingle);
			vector<int> vtemp = tess.GetNeighbors(index);
			size_t NN = vtemp.size();
			Vector2D llperiodic = 2.1 * ll - 1.1*ur;
			Vector2D urperiodic = 2.1 * ur - 1.1*ll;
			for (size_t j = 0; j < NN; ++j)
			{
				if (vtemp[j] < tess.GetPointNo())
					find_affected_cells2(tess, vtemp[j], circle, res, visited, added, periodic, toaddsingle,visited_add,ll,ur);
				else
					if (periodic)
					{
						Vector2D toadd = tess.GetMeshPoint(tess.GetOriginalIndex(vtemp[j])) -
							tess.GetMeshPoint(vtemp[j]);
						// Make sure not heading back in duplicate space
						if (!periodic_duplicate_before(added, toadd,llperiodic,urperiodic))
						{
							Vector2D oldcenter = circle.getCenter();
							circle.setCenter(toadd + oldcenter);
							toadd += toaddsingle;
							find_affected_cells2(tess, tess.GetOriginalIndex(vtemp[j]), circle, res, visited, added, periodic,
								toadd,visited_add,ll,ur);
							circle.setCenter(oldcenter);
						}
					}
			}
			return;
		}
	}
}

vector<int> find_affected_cells(const Tessellation& tess, int index, Circle& circle, vector<int> &vtemp, bool periodic,
	std::vector<Vector2D> &periodic_add)
{
	periodic_add.clear();
	vector<int> res;
	vtemp = tess.GetNeighbors(index);
	size_t N = vtemp.size();
	for (size_t i = 0; i < N; ++i)
	{

		if (vtemp[i] < tess.GetPointNo())
		{
			if (cell_circle_intersect(tess, vtemp[i], circle))
			{
				res.push_back(vtemp[i]);
				if (periodic)
					periodic_add.push_back(Vector2D());
			}
		}
		else
		{
			if (periodic)
			{
				int real_neigh = tess.GetOriginalIndex(vtemp[i]);
				Vector2D toadd = tess.GetMeshPoint(real_neigh) - tess.GetMeshPoint(vtemp[i]) ;
				circle.setCenter(circle.getCenter() + toadd);
				if (cell_circle_intersect(tess, real_neigh, circle))
				{
					res.push_back(real_neigh);
					periodic_add.push_back(toadd);
				}
				circle.setCenter(circle.getCenter() - toadd);
			}
		}
	}
	return res;
}

void find_affected_cells_recursive(const Tessellation& tess, int index, const Circle& circle,
	vector<int> &res, std::vector<Vector2D> &added, bool periodic,Vector2D const& ll,Vector2D const& ur)
{
	vector<int> visited;
	std::vector<Vector2D> visited_add;
	added.clear();
	Vector2D toaddsingle;
	find_affected_cells2(tess, index, circle, res, visited, added, periodic, toaddsingle, visited_add,ll,ur);
}
