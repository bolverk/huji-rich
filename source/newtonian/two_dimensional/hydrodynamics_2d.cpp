#include "hydrodynamics_2d.hpp"
#include "../../tessellation/calc_face_vertex_velocity.hpp"

using std::max;

namespace
{

	class EdgeLengthCalculator : public Index2Member<double>
	{
	public:

		EdgeLengthCalculator(const Tessellation& tess,
			const PhysicalGeometry& pg) :
			tess_(tess), pg_(pg) {}

		size_t getLength(void) const
		{
			return static_cast<size_t>(tess_.GetTotalSidesNumber());
		}

		double operator()(size_t i) const
		{
			return pg_.calcArea(tess_.GetEdge(static_cast<int>(i)));
		}

	private:
		const Tessellation& tess_;
		const PhysicalGeometry& pg_;
	};

	/*Vector2D AreaOverlap2(vector<Vector2D> const& poly0, vector<Vector2D> const& poly1, vector<Vector2D> const& poly0new, vector<Vector2D> const& poly1new, double R0, double R1,
	ClipperLib::Clipper &c, ClipperLib::Paths &subj, ClipperLib::Paths &clip, ClipperLib::Paths &solution)
	{
	// returns the overlap of p1 with p0. p0 with p1
	using namespace ClipperLib;
	double maxi = 0;
	for (size_t i = 0; i < poly0.size(); ++i)
	maxi = std::max(maxi, std::max(std::abs(poly0[i].x), std::abs(poly0[i].y)));
	for (size_t i = 0; i < poly1.size(); ++i)
	maxi = std::max(maxi, std::max(std::abs(poly1[i].x), std::abs(poly1[i].y)));
	int maxscale = (int)log10(maxi) + 9;

	subj[0].resize(poly0.size());
	clip[0].resize(poly1new.size());
	for (size_t i = 0; i < poly0.size(); ++i)
	//subj[0] << IntPoint((cInt)(poly0[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly0[i].y*pow(10.0, 18 - maxscale)));
	subj[0][i] = IntPoint((cInt)(poly0[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly0[i].y*pow(10.0, 18 - maxscale)));
	for (size_t i = 0; i < poly1new.size(); ++i)
	clip[0][i] = IntPoint((cInt)(poly1new[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly1new[i].y*pow(10.0, 18 - maxscale)));

	//perform intersection ...
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(clip, ptClip, true);
	c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
	Vector2D res;
	PolygonOverlap polyoverlap;
	if (!solution.empty())
	{
	if (solution[0].size() > 2)
	{
	vector<Vector2D> inter(solution[0].size());
	double factor = pow(10.0, maxscale - 18);
	for (size_t i = 0; i < solution[0].size(); ++i)
	inter[i].Set(solution[0][i].X * factor, solution[0][i].Y * factor);
	res.x = polyoverlap.PolyArea(inter);
	}
	}
	//subj[0].clear();
	//clip[0].clear();
	subj[0].resize(poly0new.size());
	clip[0].resize(poly1.size());
	c.Clear();
	solution.clear();
	for (size_t i = 0; i < poly0new.size(); ++i)
	//	subj[0] << IntPoint((cInt)(poly0new[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly0new[i].y*pow(10.0, 18 - maxscale)));
	subj[0][i] = IntPoint((cInt)(poly0new[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly0new[i].y*pow(10.0, 18 - maxscale)));
	for (size_t i = 0; i < poly1.size(); ++i)
	//clip[0] << IntPoint((cInt)(poly1[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly1[i].y*pow(10.0, 18 - maxscale)));
	clip[0][i] = IntPoint((cInt)(poly1[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly1[i].y*pow(10.0, 18 - maxscale)));
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(clip, ptClip, true);
	c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
	if (!solution.empty())
	{
	if (solution[0].size() > 2)
	{
	vector<Vector2D> inter(solution[0].size());
	double factor = pow(10.0, maxscale - 18);
	for (size_t i = 0; i < solution[0].size(); ++i)
	inter[i].Set(solution[0][i].X * factor, solution[0][i].Y * factor);
	res.y = polyoverlap.PolyArea(inter);
	}
	}
	return res;
	}

	Vector2D AreaOverlap(vector<Vector2D> const& poly0, vector<Vector2D> const& poly1, vector<Vector2D> const& poly0new, vector<Vector2D> const& poly1new, double R0, double R1)
	{
	// returns the overlap of p1 with p0. p0 with p1
	using namespace ClipperLib;
	Paths subj(1), clip(1), solution;
	double maxi = 0;
	for (size_t i = 0; i < poly0.size(); ++i)
	maxi = std::max(maxi, std::max(std::abs(poly0[i].x), std::abs(poly0[i].y)));
	for (size_t i = 0; i < poly1.size(); ++i)
	maxi = std::max(maxi, std::max(std::abs(poly1[i].x), std::abs(poly1[i].y)));
	int maxscale = (int)log10(maxi) + 9;

	subj[0].resize(poly0.size());
	clip[0].resize(poly1new.size());
	for (size_t i = 0; i < poly0.size(); ++i)
	//subj[0] << IntPoint((cInt)(poly0[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly0[i].y*pow(10.0, 18 - maxscale)));
	subj[0][i] = IntPoint((cInt)(poly0[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly0[i].y*pow(10.0, 18 - maxscale)));
	for (size_t i = 0; i < poly1new.size(); ++i)
	clip[0][i] = IntPoint((cInt)(poly1new[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly1new[i].y*pow(10.0, 18 - maxscale)));

	//perform intersection ...
	Clipper c;
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(clip, ptClip, true);
	c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
	Vector2D res;
	PolygonOverlap polyoverlap;
	if (!solution.empty())
	{
	if (solution[0].size() > 2)
	{
	vector<Vector2D> inter(solution[0].size());
	double factor = pow(10.0, maxscale - 18);
	for (size_t i = 0; i < solution[0].size(); ++i)
	inter[i].Set(solution[0][i].X * factor, solution[0][i].Y * factor);
	res.x = polyoverlap.PolyArea(inter);
	}
	}
	//subj[0].clear();
	//clip[0].clear();
	subj[0].resize(poly0new.size());
	clip[0].resize(poly1.size());
	c.Clear();
	solution.clear();
	for (size_t i = 0; i < poly0new.size(); ++i)
	//	subj[0] << IntPoint((cInt)(poly0new[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly0new[i].y*pow(10.0, 18 - maxscale)));
	subj[0][i] = IntPoint((cInt)(poly0new[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly0new[i].y*pow(10.0, 18 - maxscale)));
	for (size_t i = 0; i < poly1.size(); ++i)
	//clip[0] << IntPoint((cInt)(poly1[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly1[i].y*pow(10.0, 18 - maxscale)));
	clip[0][i] = IntPoint((cInt)(poly1[i].x*pow(10.0, 18 - maxscale)), (cInt)(poly1[i].y*pow(10.0, 18 - maxscale)));
	c.AddPaths(subj, ptSubject, true);
	c.AddPaths(clip, ptClip, true);
	c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
	if (!solution.empty())
	{
	if (solution[0].size() > 2)
	{
	vector<Vector2D> inter(solution[0].size());
	double factor = pow(10.0, maxscale - 18);
	for (size_t i = 0; i < solution[0].size(); ++i)
	inter[i].Set(solution[0][i].X * factor, solution[0][i].Y * factor);
	res.y = polyoverlap.PolyArea(inter);
	}
	}
	return res;
	}
	*/
	Vector2D PeriodDiff(Vector2D const& pointvelocity, double dx, double dy, Vector2D const& real_p,
		Vector2D const& real_p_new)
	{
		Vector2D toadd;
		if ((real_p.x - real_p_new.x)*pointvelocity.x > 0)
			toadd.x = ((real_p.x - real_p_new.x) > 0) ? dx : -dx;
		if ((real_p.y - real_p_new.y)*pointvelocity.y > 0)
			toadd.y = ((real_p.y - real_p_new.y) > 0) ? dy : -dy;
		return toadd;
	}

	/*
	Vector2D PeriodicFix1(Vector2D const& pointvelocity, double dx, double dy, Vector2D const& real_p,
		Vector2D &real_p_new, vector<Vector2D> &p0new)
	{
		Vector2D toadd;
		bool added = false;
		if ((real_p.x - real_p_new.x)*pointvelocity.x > 0)
		{
			toadd.x = ((real_p.x - real_p_new.x) > 0) ? dx : -dx;
			added = true;
		}
		if ((real_p.y - real_p_new.y)*pointvelocity.y > 0)
		{
			toadd.y = ((real_p.y - real_p_new.y) > 0) ? dy : -dy;
			added = true;
		}
		if (added)
		{
			real_p_new += toadd;
			p0new += toadd;
		}
		return toadd;
	}

	Vector2D PeriodicFix1(Vector2D const& pointvelocity, double dx, double dy, Vector2D const& real_p,
		Vector2D &real_p_new)
	{
		Vector2D toadd;
		bool added = false;
		if ((real_p.x - real_p_new.x)*pointvelocity.x > 0)
		{
			toadd.x = ((real_p.x - real_p_new.x) > 0) ? dx : -dx;
			added = true;
		}
		if ((real_p.y - real_p_new.y)*pointvelocity.y > 0)
		{
			toadd.y = ((real_p.y - real_p_new.y) > 0) ? dy : -dy;
			added = true;
		}
		if (added)
			real_p_new += toadd;
		return toadd;
	}
*/
	/*
	void PeriodicFix2(Vector2D const& real_p, Vector2D const& real_p_new, double dx, double dy, vector<Vector2D> &p1, vector<Vector2D> &p1new, Vector2D const& oldp, Vector2D const& newp,
		Vector2D &added1, Vector2D &added2)
	{
		Vector2D toadd;
		if (2 * (real_p.x - oldp.x) > dx)
			toadd += Vector2D(dx, 0);
		if (2 * (real_p.y - oldp.y) > dy)
			toadd += Vector2D(0, dy);
		if (2 * (-real_p.x + oldp.x) > dx)
			toadd += Vector2D(-dx, 0);
		if (2 * (-real_p.y + oldp.y) > dy)
			toadd += Vector2D(0, -dy);
		p1 += toadd;
		added1 = toadd;
		toadd.Set(0, 0);

		if (2 * (real_p_new.x - newp.x) > dx)
			toadd += Vector2D(dx, 0);
		if (2 * (real_p_new.y - newp.y) > dy)
			toadd += Vector2D(0, dy);
		if (2 * (-real_p_new.x + newp.x) > dx)
			toadd += Vector2D(-dx, 0);
		if (2 * (-real_p_new.y + newp.y) > dy)
			toadd += Vector2D(0, -dy);
		p1new += toadd;
		added2 = toadd;
	}
	*/
	/*
	std::pair<vector<std::pair<Vector2D, Vector2D> >, vector<std::pair<Vector2D, Vector2D> > > GetCorners(vector<vector<Vector2D> > const& oldp, vector<vector<Vector2D> > const& newp)
	{
		std::pair<vector<std::pair<Vector2D, Vector2D> >, vector<std::pair<Vector2D, Vector2D> > > res;
		res.first.resize(oldp.size());
		res.second.resize(newp.size());
		for (size_t i = 0; i < oldp.size(); ++i)
		{
			double minx(oldp[i][0].x), maxx(oldp[i][0].x), maxy(oldp[i][0].y), miny(oldp[i][0].y);
			for (size_t j = 1; j < oldp[i].size(); ++j)
			{
				minx = min(minx, oldp[i][j].x);
				maxx = max(maxx, oldp[i][j].x);
				miny = min(miny, oldp[i][j].y);
				maxy = max(maxy, oldp[i][j].y);
			}
			res.first[i].first.Set(minx, miny);
			res.first[i].second.Set(maxx, maxy);
			minx = newp[i][0].x;
			maxx = newp[i][0].x;
			maxy = newp[i][0].y;
			miny = newp[i][0].y;
			for (size_t j = 1; j < newp[i].size(); ++j)
			{
				minx = min(minx, newp[i][j].x);
				maxx = max(maxx, newp[i][j].x);
				miny = min(miny, newp[i][j].y);
				maxy = max(maxy, newp[i][j].y);
			}
			res.second[i].first.Set(minx, miny);
			res.second[i].second.Set(maxx, maxy);
		}
		return res;
	}
*/
/*
	bool NeedToCheck(std::pair<Vector2D, Vector2D> p0, std::pair<Vector2D, Vector2D> p0new, std::pair<Vector2D, Vector2D> p1, std::pair<Vector2D, Vector2D> p1new)
	{
		if (p0.first.x > p1new.second.x&&p0new.first.x > p1.second.x)
			return false;
		if (p0.second.x < p1new.first.x&&p0new.second.x < p1.first.x)
			return false;
		if (p0.first.y > p1new.second.y&&p0new.first.y > p1.second.y)
			return false;
		if (p0.second.y < p1new.first.y&&p0new.second.y < p1.first.y)
			return false;
		return true;
	}
*/
	bool DoEdgesIntersect(Edge const& e1, Edge const &e2)
	{
		TripleConstRef<Vector2D> points0(e1.vertices.first, e1.vertices.second, e2.vertices.first);
		TripleConstRef<Vector2D> points1(e1.vertices.first, e1.vertices.second, e2.vertices.second);
		if (orient2d(points0)*orient2d(points0)>0)
			return false;
		TripleConstRef<Vector2D> points2(e2.vertices.first, e2.vertices.second, e1.vertices.first);
		TripleConstRef<Vector2D> points3(e2.vertices.first, e2.vertices.second, e1.vertices.second);
		if (orient2d(points2)*orient2d(points3)>0)
			return false;
		return true;
	}

	bool AreParallel(Edge const& e1, Edge const &e2)
	{
		TripleConstRef<Vector2D> points(e2.vertices.first - e1.vertices.first,
			e2.vertices.second - e1.vertices.second, Vector2D(0, 0));
		const double small = 1e-7*e1.GetLength()*e2.GetLength();
		if (std::abs(orient2d(points))<small)
			return true;
		else
			return false;
	}

	Vector2D FindIntersection(Edge const& e1, Edge const &e2)
	{
		TripleConstRef<Vector2D> points(e1.vertices.first, e1.vertices.second, e2.vertices.first);
		const double temp = ((e2.vertices.second.x - e2.vertices.first.x)*
			(e1.vertices.second.y - e1.vertices.first.y) - (e2.vertices.second.y - e2.vertices.first.y)*
			(e1.vertices.second.x - e1.vertices.first.x));
		double alpha = orient2d(points) / temp;
		if (alpha>0.9999)
		{
			TripleConstRef<Vector2D> points2(e1.vertices.first, e1.vertices.second, e2.vertices.second);
			alpha = -orient2d(points2) / temp;
			return e2.vertices.second + alpha*(e2.vertices.first - e2.vertices.second);
		}
		return e2.vertices.first - alpha*(e2.vertices.first - e2.vertices.second);
	}

	bool RightHand(Edge const& e, Vector2D const& point)
	{
		TripleConstRef<Vector2D> points(point, e.vertices.first, e.vertices.second);
		if (orient2d(points)>0)
			return true;
		else
			return false;
	}
	
	bool IsFirstVertice(Tessellation const& tess, int cell, Vector2D const& added, Edge const& other_edge)
	{
		vector<int> const& edges = tess.GetCellEdges(cell);
		double minR_2 = tess.GetVolume(cell);
		bool res = false;
		for (size_t i = 0; i < edges.size(); ++i)
		{
			Edge const& edge_temp = tess.GetEdge(edges[i]);
			double temp = (edge_temp.vertices.first - added).distance2(other_edge.vertices.first);
			if (temp < minR_2)
			{
				minR_2 = temp;
				res = true;
			}
			temp = (edge_temp.vertices.first - added).distance2(other_edge.vertices.second);
			if (temp < minR_2)
			{
				minR_2 = temp;
				res = false;
			}
			temp = (edge_temp.vertices.second - added).distance2(other_edge.vertices.first);
			if (temp < minR_2)
			{
				minR_2 = temp;
				res = true;
			}
			temp = (edge_temp.vertices.second - added).distance2(other_edge.vertices.second);
			if (temp < minR_2)
			{
				minR_2 = temp;
				res = false;
			}
		}
		return res;
	}

	std::pair<Vector2D, Vector2D> TrianglesArea(Edge const& e1, Edge const &e2, Tessellation const& tessnew, Tessellation const& tessold,
		Vector2D const& added0, Vector2D const& added1, Vector2D const& added3, int cell_index, boost::array<Vector2D, 4> &CM)
	{ // first is the change in the old and the second the change in the new
		// x is e1.first y is e1.second
		std::pair<Vector2D, Vector2D> res;
		bool first = IsFirstVertice(tessnew, cell_index, -1 * added0, e2);
		if (e1.neighbors.first != cell_index)
			first = !first;
		Vector2D third_point;
		int sign_1;
		if (e1.neighbors.first >= 0)
		{
			third_point = first ? e2.vertices.first : e2.vertices.second;
			TripleConstRef<Vector2D> temp(e1.vertices.first, e1.vertices.second, tessold.GetMeshPoint(e1.neighbors.first));
			sign_1 = orient2d(temp) > 0 ? -1 : 1;
			TripleConstRef<Vector2D> temp2 = TripleConstRef<Vector2D>(e1.vertices.first, e1.vertices.second, third_point);
			res.first.x = 0.5*sign_1*orient2d(temp2);
			CM[0] = (e1.vertices.first + e1.vertices.second + third_point) / 3.0;
		}

		if (e1.neighbors.second >= 0)
		{
			third_point = !first ? e2.vertices.first : e2.vertices.second;
			TripleConstRef<Vector2D> temp3 = TripleConstRef<Vector2D>(e1.vertices.first, e1.vertices.second, tessold.GetMeshPoint(e1.neighbors.second));
			sign_1 = orient2d(temp3) > 0 ? -1 : 1;
			TripleConstRef<Vector2D> temp4 = TripleConstRef<Vector2D>(e1.vertices.first, e1.vertices.second, third_point);
			res.first.y = 0.5*sign_1*orient2d(temp4);
			CM[1] = (e1.vertices.first + e1.vertices.second + third_point) / 3.0;
		}

		if (e2.neighbors.first >= 0)
			first = IsFirstVertice(tessold, tessnew.GetOriginalIndex(e2.neighbors.first), added3, e1);
		else
			first = !IsFirstVertice(tessold, tessnew.GetOriginalIndex(e2.neighbors.second), added3, e1);
		third_point = first ? e1.vertices.first : e1.vertices.second;
		if (e2.neighbors.first >= 0)
		{
			TripleConstRef<Vector2D> temp5 = TripleConstRef<Vector2D>(e2.vertices.first, e2.vertices.second, tessnew.GetMeshPoint(e2.neighbors.first) - added1);
			sign_1 = orient2d(temp5) < 0 ? -1 : 1;
			TripleConstRef<Vector2D> temp6 = TripleConstRef<Vector2D>(e2.vertices.first, e2.vertices.second, third_point);
			res.second.x = 0.5*sign_1*orient2d(temp6);
			CM[2] = (e2.vertices.first + e2.vertices.second + third_point) / 3.0;
		}
		third_point = !first ? e1.vertices.first : e1.vertices.second;
		if (e2.neighbors.second >= 0)
		{
			TripleConstRef<Vector2D> temp7 = TripleConstRef<Vector2D>(e2.vertices.first, e2.vertices.second, tessnew.GetMeshPoint(e2.neighbors.second) - added1);
			sign_1 = orient2d(temp7) < 0 ? -1 : 1;
			TripleConstRef<Vector2D> temp8 = TripleConstRef<Vector2D>(e2.vertices.first, e2.vertices.second, third_point);
			res.second.y = 0.5*sign_1*orient2d(temp8);
			CM[3] = (e2.vertices.first + e2.vertices.second + third_point) / 3.0;
		}
		return res;
	}

	double EdgesArea(Edge const& e1, Edge const &e2, Vector2D const& pointold,
		Vector2D const& pointnew)
	{
		// e1 is old tess, e2 is new tess
		if (!DoEdgesIntersect(e1, e2))
		{
			double res = 0;
			if (!RightHand(e1, pointold))
			{
				if (!RightHand(e2, pointnew))
				{
					TripleConstRef<Vector2D> temp(e1.vertices.second, e2.vertices.second, e2.vertices.first);
					res += orient2d(temp);
					TripleConstRef<Vector2D> temp2(e1.vertices.second, e2.vertices.first, e1.vertices.first);
					res += orient2d(temp2);
				}
				else
				{
					TripleConstRef<Vector2D> temp(e1.vertices.second, e2.vertices.first, e2.vertices.second);
					res += orient2d(temp);
					TripleConstRef<Vector2D> temp2(e1.vertices.second, e2.vertices.second, e1.vertices.first);
					res += orient2d(temp2);
				}
			}
			else
			{
				if (!RightHand(e2, pointnew))
				{
					TripleConstRef<Vector2D> temp(e1.vertices.first, e2.vertices.second, e2.vertices.first);
					res += orient2d(temp);
					TripleConstRef<Vector2D> temp2(e1.vertices.first, e2.vertices.first, e1.vertices.second);
					res += orient2d(temp2);
				}
				else
				{
					TripleConstRef<Vector2D> temp(e1.vertices.first, e2.vertices.first, e2.vertices.second);
					res += orient2d(temp);
					TripleConstRef<Vector2D> temp2(e1.vertices.first, e2.vertices.second, e1.vertices.second);
					res += orient2d(temp2);
				}
			}
			return 0.5*res;
		}
		else
		{
			if (AreParallel(e1, e2))
				return 0;
			Vector2D inter = FindIntersection(e1, e2);
			double res = 0;
			if (!RightHand(e1, pointold))
			{
				if (!RightHand(e2, pointnew))
				{
					TripleConstRef<Vector2D> temp(inter, e1.vertices.second, e2.vertices.second);
					res += orient2d(temp);
					TripleConstRef<Vector2D> temp2(inter, e2.vertices.first, e1.vertices.first);
					res += orient2d(temp2);
				}
				else
				{
					TripleConstRef<Vector2D> temp(inter, e1.vertices.second, e2.vertices.first);
					res += orient2d(temp);
					TripleConstRef<Vector2D> temp2(inter, e2.vertices.second, e1.vertices.first);
					res += orient2d(temp2);
				}
			}
			else
			{
				if (!RightHand(e2, pointnew))
				{
					TripleConstRef<Vector2D> temp(inter, e1.vertices.first, e2.vertices.first);
					res += orient2d(temp);
					TripleConstRef<Vector2D> temp2(inter, e2.vertices.second, e1.vertices.second);
					res += orient2d(temp2);
				}
				else
				{
					TripleConstRef<Vector2D> temp(inter, e1.vertices.first, e2.vertices.second);
					res += orient2d(temp);
					TripleConstRef<Vector2D> temp2(inter, e2.vertices.first, e1.vertices.second);
					res += orient2d(temp2);
				}
			}
			return res*0.5;
		}
	}
	/*
	int GetEdgeIndex(size_t i, Tessellation const& tessold, Tessellation const& tessnew,
		bool first = true)
	{
		Edge const& eold = tessold.GetEdge((int)i);
		bool usefirst = (eold.neighbors.first >= 0 && eold.neighbors.first<tessold.GetPointNo())
			? true : false;
		int n0 = tessold.GetOriginalIndex(eold.neighbors.first);
		int n1 = tessold.GetOriginalIndex(eold.neighbors.second);
		if ((usefirst&&first) || (!first&&!usefirst))
		{
			vector<int> const& edges = tessnew.GetCellEdges(n0);
			for (size_t j = 0; j<edges.size(); ++j)
			{
				Edge const& etemp = tessnew.GetEdge(edges[j]);
				if ((tessnew.GetOriginalIndex(etemp.neighbors.first) == n0&&
					tessnew.GetOriginalIndex(etemp.neighbors.second) == n1) ||
					(tessnew.GetOriginalIndex(etemp.neighbors.first) == n1&&
					tessnew.GetOriginalIndex(etemp.neighbors.second) == n0))
					return edges[j];
			}
		}
		else
		{
			vector<int> const& edges = tessnew.GetCellEdges(n1);
			for (size_t j = 0; j<edges.size(); ++j)
			{
				Edge const& etemp = tessnew.GetEdge(edges[j]);
				if ((tessnew.GetOriginalIndex(etemp.neighbors.first) == n0&&
					tessnew.GetOriginalIndex(etemp.neighbors.second) == n1) ||
					(tessnew.GetOriginalIndex(etemp.neighbors.first) == n1&&
					tessnew.GetOriginalIndex(etemp.neighbors.second) == n0))
					return edges[i];
			}
		}
		return -1;
	}
*/
	int GetEdgeIndex(Tessellation const& tessnew, int n0, int n1, int cell_index)
	{
		if (cell_index < 0)
			return -1;
		vector<int> const& edges = tessnew.GetCellEdges(cell_index);
		for (size_t j = 0; j<edges.size(); ++j)
		{
			Edge const& etemp = tessnew.GetEdge(edges[j]);
			if ((tessnew.GetOriginalIndex(etemp.neighbors.first) == n0&&
				tessnew.GetOriginalIndex(etemp.neighbors.second) == n1) ||
				(tessnew.GetOriginalIndex(etemp.neighbors.first) == n1&&
				tessnew.GetOriginalIndex(etemp.neighbors.second) == n0))
				return edges[j];
		}
		return -1;
	}

	Primitive GetDonorPrimitive(Primitive const& c0, Primitive const& c1,
		bool seconddonor, EquationOfState const& eos)
	{
		Primitive p_mid;
		double ratio = c1.Density / c0.Density;
		if (ratio<2 && ratio>0.5)
			p_mid.Density = 0.5*(c1.Density + c0.Density);
		else
		{
			if (seconddonor)
				p_mid.Density = c1.Density;
			else
				p_mid.Density = c0.Density;
		}
		ratio = c1.Pressure / c0.Pressure;
		if (ratio<2 && ratio>0.5)
			p_mid.Pressure = 0.5*(c1.Pressure + c0.Pressure);
		else
		{
			if (seconddonor)
				p_mid.Pressure = c1.Pressure;
			else
				p_mid.Pressure = c0.Pressure;
		}
		ratio = c1.Velocity.x / c0.Velocity.x;
		if (ratio<2 && ratio>0.5)
			p_mid.Velocity.x = 0.5*(c1.Velocity.x + c0.Velocity.x);
		else
		{
			if (seconddonor)
				p_mid.Velocity.x = c1.Velocity.x;
			else
				p_mid.Velocity.x = c0.Velocity.x;
		}
		ratio = c1.Velocity.y / c0.Velocity.y;
		if (ratio<2 && ratio>0.5)
			p_mid.Velocity.y = 0.5*(c1.Velocity.y + c0.Velocity.y);
		else
		{
			if (seconddonor)
				p_mid.Velocity.y = c1.Velocity.y;
			else
				p_mid.Velocity.y = c0.Velocity.y;
		}
		p_mid.Energy = eos.dp2e(p_mid.Density, p_mid.Pressure);
		return p_mid;
	}


	bool SameNeighbors(Edge const& e1, Edge const& e2, Tessellation const& tess)
	{
		if ((tess.GetOriginalIndex(e1.neighbors.first) == tess.GetOriginalIndex(e2.neighbors.first))
			&& (tess.GetOriginalIndex(e1.neighbors.second) ==
			tess.GetOriginalIndex(e2.neighbors.second)))
			return true;
		if ((tess.GetOriginalIndex(e1.neighbors.second) == tess.GetOriginalIndex(e2.neighbors.first))
			&& (tess.GetOriginalIndex(e1.neighbors.first) ==
			tess.GetOriginalIndex(e2.neighbors.second)))
			return true;
		return false;
	}

	void GetAdjacentNeigh(Tessellation const& tessold, Edge const& edge, int cell_index,
		Vector2D const& added, int &n0other, int &n1other)
	{
		vector<int> const& edges = tessold.GetCellEdges(cell_index);
		vector<double> dist;
		Vector2D point;
		dist.reserve(edges.size());
		for (size_t i = 0; i<edges.size(); ++i)
		{
			Edge const& e = tessold.GetEdge(edges[i]);
			if (SameNeighbors(edge, e, tessold))
				dist.push_back(edge.GetLength() * 10);
			else
				dist.push_back(min(min(min(abs(e.vertices.first - edge.vertices.first + added),
				abs(e.vertices.second - edge.vertices.first + added)),
				abs(e.vertices.first - edge.vertices.second + added)),
				abs(e.vertices.second - edge.vertices.second + added)));
		}
		vector<size_t> indeces = sort_index(dist);
		// The location of the edges to the left and the right of the old edge
		n0other = (tessold.GetOriginalIndex(tessold.GetEdge(edges[indeces[0]]).neighbors.first) == cell_index) ?
			tessold.GetOriginalIndex(tessold.GetEdge(edges[indeces[0]]).neighbors.second) :
			tessold.GetOriginalIndex(tessold.GetEdge(edges[indeces[0]]).neighbors.first);
		n1other = (tessold.GetOriginalIndex(tessold.GetEdge(edges[indeces[1]]).neighbors.first) == cell_index) ?
			tessold.GetOriginalIndex(tessold.GetEdge(edges[indeces[1]]).neighbors.second) :
			tessold.GetOriginalIndex(tessold.GetEdge(edges[indeces[1]]).neighbors.first);
	}

	// Find the new edge in the edge flip
	int NewEdgeIndex(Tessellation const& tessold, Tessellation const& tessnew, int cell_index, Edge const& edge, int other_index)
	{
		int n0, n1;
		int res = -1;
		// Find n0 and n1 the new flip neighbors
		GetAdjacentNeigh(tessold, edge, cell_index, Vector2D(0, 0), n0, n1);
		if (n0 == -1)
		{
			n0 = n1;
			n1 = -1;
		}
		vector<int> const& edges = tessnew.GetCellEdges(n0);
		for (size_t i = 0; i < edges.size(); ++i)
		{
			Edge const& edge_temp = tessnew.GetEdge(edges[i]);
			int n0new = tessnew.GetOriginalIndex(edge_temp.neighbors.first), n1new = tessnew.GetOriginalIndex(edge_temp.neighbors.second);
			if (n0new == n1 || n1new == n1)
				res = edges[i];
		}
		if (res == -1)
			return -1;
		int n0new, n1new;
		GetAdjacentNeigh(tessnew, tessnew.GetEdge(res), n0, Vector2D(0, 0), n0new, n1new);
		if ((n0new == cell_index&&n1new == other_index) || (n1new == cell_index&&other_index == n0new))
			return res;
		// We didn't find a new edge
		return -1;
	}
}

vector<Conserved> FluxFix2(Tessellation const& tessold, Tessellation const& tessmid,
	Tessellation const& tessnew, vector<Vector2D> const& pointvelocity, double dt,
	vector<Primitive> const& cells, vector<Conserved> &/*fluxes*/,
	vector<Vector2D> const& fv, OuterBoundary const& outer, EquationOfState const& eos,
	HydroBoundaryConditions const& hbc, vector<vector<double> > const& tracers, vector<vector<double> > &dtracers)
{
	/*voronoi_loggers::BinLogger log("c:\\v0.bin"), log2("c:\\v1.bin"), log3("c:\\v2.bin");
	log.output(tessold);
	log2.output(tessnew);
	log3.output(tessmid);*/
	bool traceractive = !tracers.empty();
	size_t tracer_dim = 0;
	dtracers.clear();
	if (traceractive)
	{
		tracer_dim = tracers[0].size();
		vector<double> ttemp(tracer_dim, 0);
		dtracers.resize(tracers.size(), ttemp);
	}
	int npoints = tessold.GetPointNo();
	vector<Conserved> res(static_cast<size_t>(npoints));
	vector<vector<Vector2D> > chull_old(static_cast<size_t>(npoints)),
	  chull_new(static_cast<size_t>(npoints));
	const double dx = outer.GetGridBoundary(Right) - outer.GetGridBoundary(Left);
	const double dy = outer.GetGridBoundary(Up) - outer.GetGridBoundary(Down);
	// Fix the fluxes
	Vector2D temp0, temp1;
	std::set<std::pair<int, int > > flipped_set, only_mid;
	size_t nedgesold = static_cast<size_t>(tessold.GetTotalSidesNumber());
	for (size_t i = 0; i < nedgesold; ++i)
	{
	  Edge edge = tessold.GetEdge(static_cast<int>(i));

		bool rigid_edge = false;
		if (edge.neighbors.first < 0 || edge.neighbors.second < 0)
			rigid_edge = true;

		int cell_index, other_index;
		if (hbc.IsGhostCell(edge.neighbors.first, tessold))
		{
			cell_index = edge.neighbors.second;
			other_index = tessold.GetOriginalIndex(edge.neighbors.first);
		}
		else
		{
			cell_index = edge.neighbors.first;
			other_index = tessold.GetOriginalIndex(edge.neighbors.second);
		}
		int mid_index = GetEdgeIndex(tessmid, cell_index, other_index, cell_index);
		int new_index = GetEdgeIndex(tessnew, cell_index, other_index, cell_index);
		double dA_flux = 0;

		if (mid_index >= 0 && !rigid_edge)
		{
			Edge const& edge2 = tessmid.GetEdge(mid_index);
			const Vector2D norm = tessmid.GetMeshPoint(edge2.neighbors.second) - tessmid.GetMeshPoint(edge2.neighbors.first);
			dA_flux = ScalarProd(norm, fv[static_cast<size_t>(mid_index)])*dt*edge2.GetLength() / abs(norm);
			if (tessmid.GetOriginalIndex(edge2.neighbors.first) == other_index)
				dA_flux *= -1;
		}

		Vector2D real_p, real_p_new;
		real_p = tessold.GetMeshPoint(cell_index);
		real_p_new = tessnew.GetMeshPoint(cell_index);
		Vector2D added(0, 0);
		// check if periodic jumped
		if (outer.GetBoundaryType() == Periodic)
		  added = PeriodDiff(pointvelocity[static_cast<size_t>(cell_index)], dx, dy, real_p, real_p_new);

		// Did the edge disappear?
		if (new_index<0)
		{
			// Is there a corresponding new edge?
			int new_edge = NewEdgeIndex(tessold, tessnew, cell_index, edge, other_index);
			if (new_edge < 0)
				continue;

			boost::array<Vector2D, 4> CM;
			Edge edge_new = tessnew.GetEdge(new_edge);
			// Did it jump?
			Vector2D addedother;
			vector<int> const& cell_edges = tessold.GetCellEdges(cell_index);
			if (outer.GetBoundaryType() == Periodic)
			{
				for (size_t j = 0; j < cell_edges.size(); ++j)
				{
					Edge const& edge_temp = tessold.GetEdge(cell_edges[j]);
					if (tessold.GetOriginalIndex(edge_temp.neighbors.first) == tessnew.GetOriginalIndex(edge_new.neighbors.first))

					{
						addedother = tessold.GetMeshPoint(tessnew.GetOriginalIndex(edge_new.neighbors.first)) - tessold.GetMeshPoint(edge_temp.neighbors.first);
						break;
					}
					if (tessold.GetOriginalIndex(edge_temp.neighbors.second) == tessnew.GetOriginalIndex(edge_new.neighbors.first))
					{
						addedother = tessold.GetMeshPoint(tessnew.GetOriginalIndex(edge_new.neighbors.first)) - tessold.GetMeshPoint(edge_temp.neighbors.second);
						break;
					}
				}
			}
			Vector2D added_edge;
			if (outer.GetBoundaryType() == Periodic)
			{
				vector<int> const& new_edges = tessnew.GetCellEdges(cell_index);
				for (size_t j = 0; j < new_edges.size(); ++j)
				{
					Edge const& edge_temp = tessnew.GetEdge(new_edges[j]);
					if (tessnew.GetOriginalIndex(edge_temp.neighbors.first) == tessnew.GetOriginalIndex(edge_new.neighbors.first))
					{
						added_edge = tessnew.GetMeshPoint(edge_new.neighbors.first) - tessnew.GetMeshPoint(edge_temp.neighbors.first);
						break;
					}
					if (tessnew.GetOriginalIndex(edge_temp.neighbors.second) == tessnew.GetOriginalIndex(edge_new.neighbors.first))
					{
						added_edge = tessnew.GetMeshPoint(edge_new.neighbors.first) - tessnew.GetMeshPoint(edge_temp.neighbors.second);
						break;
					}
				}
			}
			added_edge -= added;
			edge_new.vertices.first -= added_edge;
			edge_new.vertices.second -= added_edge;
			std::pair<Vector2D, Vector2D> areas = TrianglesArea(edge, edge_new, tessnew, tessold, added, added_edge, addedother, cell_index, CM);
			vector<double> dA(4);
			if (mid_index >= 0 && tessmid.GetOriginalIndex(tessmid.GetEdge(mid_index).neighbors.first) == other_index)
				dA_flux *= -1;
			if (mid_index >= 0 && tessmid.GetOriginalIndex(tessmid.GetEdge(mid_index).neighbors.first) == tessold.GetOriginalIndex(edge.neighbors.second))
				dA_flux *= -1;
			dA[0] = areas.first.x - dA_flux;
			dA[1] = areas.first.y + dA_flux;
			dA[2] = areas.second.x;
			dA[3] = areas.second.y;


			int mid_index2 = GetEdgeIndex(tessmid, tessnew.GetOriginalIndex(edge_new.neighbors.first),
				tessnew.GetOriginalIndex(edge_new.neighbors.second), tessnew.GetOriginalIndex(edge_new.neighbors.first));

			if (edge_new.neighbors.first == -1 || edge_new.neighbors.second == -1)
				mid_index2 = -1;

			double dA_flux2 = 0;

			if (mid_index2 >= 0)
			{
				Edge const& edge2 = tessmid.GetEdge(mid_index2);
				const Vector2D norm = tessmid.GetMeshPoint(edge2.neighbors.second) - tessmid.GetMeshPoint(edge2.neighbors.first);
				dA_flux2 = ScalarProd(norm, fv[static_cast<size_t>(mid_index2)])*dt*edge2.GetLength() / abs(norm);
				if (tessmid.GetOriginalIndex(edge2.neighbors.first) != tessnew.GetOriginalIndex(edge_new.neighbors.first))
					dA_flux2 *= -1;
				dA[2] -= dA_flux2;
				dA[3] += dA_flux2;
			}

			// The conserved to remove
			Conserved TotalRemoved;
			vector<double> TracerRemoved(tracer_dim, 0);
			vector<double> ttremove;
			double total_added_area = 0;
			if (dA[0] < 0)
			{
				if (traceractive)
				  ttremove = -dA[0] * cells[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.first))].Density*
				    tracers[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.first))];
				Conserved toremove = -dA[0] * Primitive2Conserved(cells[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.first))]);
				if (!hbc.IsGhostCell(edge.neighbors.first, tessold))
				{
				  res[static_cast<size_t>(edge.neighbors.first)] -= toremove;
					if (traceractive)
					  dtracers[static_cast<size_t>(edge.neighbors.first)] = dtracers[static_cast<size_t>(edge.neighbors.first)] - ttremove;
				}
				TotalRemoved += toremove;
				if (traceractive)
					TracerRemoved = TracerRemoved + ttremove;
			}
			else
				total_added_area += dA[0];
			if (dA[1] < 0)
			{
				if (traceractive)
				  ttremove = -dA[1] * cells[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.second))].Density*
				    tracers[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.second))];
				Conserved toremove = -dA[1] * Primitive2Conserved(cells[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.second))]);
				if (!hbc.IsGhostCell(edge.neighbors.second, tessold))
				{
				  res[static_cast<size_t>(edge.neighbors.second)] -= toremove;
					if (traceractive)
					  dtracers[static_cast<size_t>(edge.neighbors.second)] = dtracers[static_cast<size_t>(edge.neighbors.second)] - ttremove;
				}
				TotalRemoved += toremove;
				if (traceractive)
					TracerRemoved = TracerRemoved + ttremove;
			}
			else
				total_added_area += dA[1];

			bool both_real = !hbc.IsGhostCell(edge.neighbors.first, tessold) && !hbc.IsGhostCell(edge.neighbors.second, tessold);
			std::pair<int, int> neighbors(tessold.GetOriginalIndex(edge.neighbors.first), tessold.GetOriginalIndex(edge.neighbors.second));
			bool need_to_calc = both_real || (flipped_set.count(neighbors) == 0);

			if (dA[2] < 0)
			{
				if (traceractive)
				  ttremove = -dA[2] * cells[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))].Density*
				    tracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))];
				Conserved toremove = -dA[2] * Primitive2Conserved(cells[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))]);
				if (need_to_calc&&tessnew.GetOriginalIndex(edge_new.neighbors.first) >= 0)
				{
				  res[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))] -= toremove;
					if (traceractive)
					  dtracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))] = dtracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))] -
						ttremove;
				}
				TotalRemoved += toremove;
				if (traceractive)
					TracerRemoved = TracerRemoved + ttremove;
			}
			else
				total_added_area += dA[2];
			if (dA[3] < 0)
			{
				if (traceractive)
				  ttremove = -dA[3] * cells[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))].Density*
				    tracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))];
				Conserved toremove = -dA[3] * Primitive2Conserved(cells[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))]);
				if (need_to_calc&&tessnew.GetOriginalIndex(edge_new.neighbors.second) >= 0)
				{
				  res[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))] -= toremove;
					if (traceractive)
					  dtracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))] = dtracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))] -
						ttremove;
				}
				TotalRemoved += toremove;
				if (traceractive)
					TracerRemoved = TracerRemoved + ttremove;
			}
			else
				total_added_area += dA[3];

			// Add the conserved
			if (dA[0] > 0)
				if (!hbc.IsGhostCell(edge.neighbors.first, tessold))
				{
				  res[static_cast<size_t>(edge.neighbors.first)] += (dA[0] / total_added_area)*TotalRemoved;
					if (traceractive)
					  dtracers[static_cast<size_t>(edge.neighbors.first)] = dtracers[static_cast<size_t>(edge.neighbors.first)] + (dA[0] / total_added_area)*TracerRemoved;
				}
			if (dA[1] > 0)
				if (!hbc.IsGhostCell(edge.neighbors.second, tessold))
				{
				  res[static_cast<size_t>(edge.neighbors.second)] += (dA[1] / total_added_area)*TotalRemoved;
					if (traceractive)
					  dtracers[static_cast<size_t>(edge.neighbors.second)] = dtracers[static_cast<size_t>(edge.neighbors.second)] + (dA[1] / total_added_area)*TracerRemoved;
				}
			if (dA[2] > 0)
				if (need_to_calc&&tessnew.GetOriginalIndex(edge_new.neighbors.first) >= 0)
				{
				  res[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))] += (dA[2] / total_added_area)*TotalRemoved;
					if (traceractive)
					  dtracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))] = dtracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))]
						+ (dA[2] / total_added_area)*TracerRemoved;
				}
			if (dA[3] > 0)
				if (need_to_calc&&tessnew.GetOriginalIndex(edge_new.neighbors.second) >= 0)
				{
				  res[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))] += (dA[3] / total_added_area)*TotalRemoved;
					if (traceractive)
					  dtracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))] = dtracers[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))]
						+ (dA[3] / total_added_area)*TracerRemoved;
				}

			if (need_to_calc && !both_real)
			{
				flipped_set.insert(neighbors);
				std::pair<int, int> neigh_temp(neighbors.second, neighbors.first);
				flipped_set.insert(neigh_temp);
			}
			continue;
		}
		if (other_index == -1)
			continue;
		Edge edge_other = tessnew.GetEdge(new_index);
		if (outer.GetBoundaryType() == Periodic)
		{
			edge_other.vertices.first += added;
			edge_other.vertices.second += added;
			real_p_new += added;
		}
		double area = EdgesArea(edge, edge_other, real_p, real_p_new);
		Primitive p_mid = GetDonorPrimitive(cells[static_cast<size_t>(cell_index)], cells[static_cast<size_t>(other_index)],
			(area - dA_flux)>0, eos);
		//Primitive p_mid = (-dA_flux + area) > 0 ? cells[other_index] : cells[cell_index];
		Conserved toadd((-dA_flux + area)*Primitive2Conserved(p_mid));
		vector<double> ttoadd;
		if (traceractive)
		{
		  ttoadd = (-dA_flux + area) > 0 ? tracers[static_cast<size_t>(other_index)] : tracers[static_cast<size_t>(cell_index)];
			ttoadd = (-dA_flux + area)*p_mid.Density*ttoadd;
		}
		if (other_index == edge.neighbors.first || other_index == edge.neighbors.second)
		{
		  res[static_cast<size_t>(other_index)] -= toadd;
			if (traceractive)
			{
			  dtracers[static_cast<size_t>(other_index)] = dtracers[static_cast<size_t>(other_index)] - ttoadd;
			}
		}
		if (cell_index == edge.neighbors.first || cell_index == edge.neighbors.second)
		{
		  res[static_cast<size_t>(cell_index)] += toadd;
			if (traceractive)
			{
			  dtracers[static_cast<size_t>(cell_index)] = dtracers[static_cast<size_t>(cell_index)] + ttoadd;
			}
		}

		if (mid_index < 0)
		{
			// Was there a flip in mid step?
			int mid_edge = NewEdgeIndex(tessold, tessmid, cell_index, edge, other_index);
			if (mid_edge < 0)
				continue;
			Edge const& edge2 = tessmid.GetEdge(mid_edge);
			if (edge2.neighbors.first < 0 || edge2.neighbors.second < 0)
				continue;
			const Vector2D norm = tessmid.GetMeshPoint(edge2.neighbors.second) - tessmid.GetMeshPoint(edge2.neighbors.first);
			double dA_flux2 = ScalarProd(norm, fv[static_cast<size_t>(mid_edge)])*dt*edge2.GetLength() / abs(norm);
			Conserved toadd2;
			if (dA_flux2 > 0)
			{
			  toadd2 = dA_flux2*Primitive2Conserved(cells[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.second))]);
				if (traceractive)
				  ttoadd = dA_flux2*cells[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.second))].Density*
				    tracers[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.second))];
			}
			else
			{
			  toadd2 = dA_flux2*Primitive2Conserved(cells[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.first))]);
				if (traceractive)
				  ttoadd = dA_flux2*cells[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.first))].Density*
				    tracers[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.first))];
			}
			std::pair<int, int> mid_neigh(tessmid.GetOriginalIndex(edge2.neighbors.first), tessmid.GetOriginalIndex(edge2.neighbors.second));
			std::pair<int, int> mid_neigh2(mid_neigh.second, mid_neigh.first);
			if (only_mid.count(mid_neigh) == 0)
			{
				only_mid.insert(mid_neigh);
				only_mid.insert(mid_neigh2);
				res[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.second))] += toadd2;
				res[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.first))] -= toadd2;
				if (traceractive)
				{
				  dtracers[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.second))] = dtracers[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.second))]
						+ ttoadd;
				  dtracers[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.first))] = dtracers[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.first))]
						- ttoadd;
				}
			}
		}
	}
	return res;
}


/*
vector<Conserved> FluxFix(Tessellation const& tessold, Tessellation const& tessmid,
Tessellation const& tessnew, vector<Vector2D> const& pointvelocity, double dt,
vector<Primitive> const& cells, vector<Conserved> const& fluxes,
vector<Vector2D> const& fv, OuterBoundary const& outer, EquationOfState const& eos,
HydroBoundaryConditions const& hbc)
{
int npoints = tessold.GetPointNo();
vector<Conserved> res((size_t)npoints);
vector<vector<Vector2D> > chull_old(npoints), chull_new(npoints);
const double dx = outer.GetGridBoundary(Right) - outer.GetGridBoundary(Left);
const double dy = outer.GetGridBoundary(Up) - outer.GetGridBoundary(Down);
// Create convex hull data
for (int i = 0; i < npoints; ++i)
{
ConvexHull(chull_old[(size_t)i], &tessold, i);
ConvexHull(chull_new[(size_t)i], &tessnew, i);
}
std::pair<vector<std::pair<Vector2D, Vector2D> >, vector<std::pair<Vector2D, Vector2D> > > corners = GetCorners(chull_old, chull_new);
// Now fix the neighbors
Vector2D periodicadded, periodicadded1, periodicadded2;
for (int i = 0; i < npoints; ++i)
{
vector<int> nneigh, neigh, neigh_new;
neigh = tessmid.GetNeighbors(i);
for (size_t j = 0; j < neigh.size(); ++j)
neigh[j] = tessmid.GetOriginalIndex(neigh[j]);
sort(neigh.begin(), neigh.end());
neigh_new = tessnew.GetNeighbors(i);
for (size_t j = 0; j < neigh_new.size(); ++j)
neigh_new[j] = tessnew.GetOriginalIndex(neigh_new[j]);
sort(neigh_new.begin(), neigh_new.end());
// Did we change neigbhors?
//		if(neigh_new==neigh)
//			continue;
tessmid.GetNeighborNeighbors(nneigh, i);
vector<Vector2D> p0(chull_old[i]), p0new(chull_new[i]);
Vector2D real_p = tessold.GetMeshPoint(i);
Vector2D real_p_new = tessnew.GetMeshPoint(i);
const double R0 = tessold.GetWidth(i);
// check if periodic jumped
if (outer.GetBoundaryType() == Periodic)
periodicadded = PeriodicFix1(pointvelocity[i], dx, dy, real_p, real_p_new, p0new);
// Can do this better by checking if neighbor in new is not in old neighbor list
for (size_t j = 0; j < nneigh.size(); ++j)
{
if (nneigh[j] < 0)
continue;
if (!std::binary_search(neigh.begin(), neigh.end(), nneigh[j]))
{
vector<Vector2D> p1(chull_old[nneigh[j]]), p1new(chull_new[nneigh[j]]);
// check if periodic
if (outer.GetBoundaryType() == Periodic)
PeriodicFix2(real_p, real_p_new, dx, dy, p1, p1new, tessold.GetMeshPoint(nneigh[j]), tessnew.GetMeshPoint(nneigh[j]), periodicadded1, periodicadded2);
if (!NeedToCheck(corners.first[i], std::pair<Vector2D, Vector2D>(corners.second[i].first + periodicadded, corners.second[i].second + periodicadded),
std::pair<Vector2D, Vector2D>(corners.first[nneigh[j]].first + periodicadded1, corners.first[nneigh[j]].second + periodicadded1),
std::pair<Vector2D, Vector2D>(corners.second[nneigh[j]].first + periodicadded2, corners.second[nneigh[j]].second + periodicadded2)))
continue;
Vector2D dA = AreaOverlap(p0, p1, p0new, p1new, R0, tessold.GetWidth(nneigh[j]));
res[i] += dA.y*Primitive2Conserved(cells[nneigh[j]]);
res[i] -= dA.x*Primitive2Conserved(cells[i]);
}
}
}

ClipperLib::Clipper c;
ClipperLib::Paths subj(1), clip(1), solution;

// Fix the fluxes
Vector2D temp0, temp1;
for (size_t i = 0; i < fluxes.size(); ++i)
{
Edge const& edge = tessmid.GetEdge((int)i);
if (edge.neighbors.first < 0 || edge.neighbors.second < 0)
continue;
int n0 = tessmid.GetOriginalIndex(edge.neighbors.first);
int n1 = tessmid.GetOriginalIndex(edge.neighbors.second);
const Vector2D norm = tessmid.GetMeshPoint(edge.neighbors.second) - tessmid.GetMeshPoint(edge.neighbors.first);
const double dA_flux = ScalarProd(norm, fv[i])*dt*edge.GetLength() / abs(norm);
Vector2D real_p, real_p_new;
real_p = tessold.GetMeshPoint(n0);
real_p_new = tessnew.GetMeshPoint(n0);
vector<Vector2D> p0(chull_old[n0]), p0new(chull_new[n0]);
vector<Vector2D> p1(chull_old[n1]), p1new(chull_new[n1]);
// check if periodic jumped
if (outer.GetBoundaryType() == Periodic)
{
PeriodicFix1(pointvelocity[n0], dx, dy, real_p, real_p_new, p0new);
PeriodicFix2(real_p, real_p_new, dx, dy, p1, p1new, tessold.GetMeshPoint(n1), tessnew.GetMeshPoint(n1), temp0, temp1);
}
const double R0 = tessold.GetWidth(n0);
c.Clear();
Vector2D dA = AreaOverlap2(p0, p1, p0new, p1new, R0, tessold.GetWidth((n1)),c,subj,clip,solution);

Primitive p_mid;
double ratio = cells[n1].Density / cells[n0].Density;
if (ratio<2 && ratio>0.5)
p_mid.Density = 0.5*(cells[n1].Density + cells[n0].Density);
else
{
if ((dA_flux - dA.y + dA.x)>0)
p_mid.Density = cells[n1].Density;
else
p_mid.Density = cells[n0].Density;
}
ratio = cells[n1].Pressure / cells[n0].Pressure;
if (ratio<2 && ratio>0.5)
p_mid.Pressure = 0.5*(cells[n1].Pressure + cells[n0].Pressure);
else
{
if ((dA_flux - dA.y + dA.x)>0)
p_mid.Pressure = cells[n1].Pressure;
else
p_mid.Pressure = cells[n0].Pressure;
}
ratio = cells[n1].Velocity.x / cells[n0].Velocity.x;
if (ratio<2 && ratio>0.5)
p_mid.Velocity.x = 0.5*(cells[n1].Velocity.x + cells[n0].Velocity.x);
else
{
if ((dA_flux - dA.y + dA.x)>0)
p_mid.Velocity.x = cells[n1].Velocity.x;
else
p_mid.Velocity.x = cells[n0].Velocity.x;
}
ratio = cells[n1].Velocity.y / cells[n0].Velocity.y;
if (ratio<2 && ratio>0.5)
p_mid.Velocity.y = 0.5*(cells[n1].Velocity.y + cells[n0].Velocity.y);
else
{
if ((dA_flux - dA.y + dA.x)>0)
p_mid.Velocity.y = cells[n1].Velocity.y;
else
p_mid.Velocity.y = cells[n0].Velocity.y;
}
p_mid.Energy = eos.dp2e(p_mid.Density, p_mid.Pressure);
Conserved toadd((dA_flux - dA.y + dA.x)*Primitive2Conserved(p_mid));
if (!hbc.IsGhostCell(edge.neighbors.second, tessmid))
res[n1] += toadd;
if (!hbc.IsGhostCell(edge.neighbors.first, tessmid))
res[n0] -= toadd;
}
return res;
}
*/

FaceVertexVelocityCalculator::FaceVertexVelocityCalculator
(const Tessellation& tess,
const vector<Vector2D>& point_velocities,
const Vector2D std::pair<Vector2D, Vector2D>::* const member,
const vector<Vector2D>& control,
const HydroBoundaryConditions& hbc) :
tess_(tess),
point_velocities_(point_velocities),
member_(member),
control_(control),
hbc_(hbc) {}

size_t FaceVertexVelocityCalculator::getLength(void) const
{
	return static_cast<size_t>(tess_.GetTotalSidesNumber());
}

Vector2D FaceVertexVelocityCalculator::operator()(size_t i) const
{
	const Edge& edge = tess_.GetEdge(static_cast<int>(i));
	if (hbc_.IsBoundary(edge, tess_))
		return control_[i];
	return calc_face_vertex_velocity
		(tess_.GetMeshPoint(edge.neighbors.first),
		point_velocities_[static_cast<size_t>(edge.neighbors.first)],
		tess_.GetMeshPoint(edge.neighbors.second),
		point_velocities_[static_cast<size_t>(edge.neighbors.second)],
		edge.vertices.*member_);
}

int get_other_index(const Edge& edge, const int index)
{
	if (edge.neighbors.first == index && edge.neighbors.second != index)
		return edge.neighbors.second;
	else if (edge.neighbors.second == index && edge.neighbors.first != index)
		return edge.neighbors.first;
	else
		throw UniversalError("Something went wrong in Hydrodynamics2D::get_other_index");
}

namespace {
	Vector2D calc_representing_point(Tessellation const& tess,
		int index,
		bool cm_flag)
	{
		if (cm_flag)
			return tess.GetCellCM(index);
		else
			return tess.GetMeshPoint(index);
	}

	Primitive initialize_single_cell(Tessellation const& tess,
		int index,
		bool cm_flag,
		SpatialDistribution const& density,
		SpatialDistribution const& pressure,
		SpatialDistribution const& xvelocity,
		SpatialDistribution const& yvelocity,
		EquationOfState const& eos)
	{
		const Vector2D r = calc_representing_point(tess, index, cm_flag);
		return CalcPrimitive(density(r),
			pressure(r),
			Vector2D(xvelocity(r),
			yvelocity(r)),
			eos);
	}

	class CellInitializer : public Index2Member<Primitive>
	{
	public:

		CellInitializer(Tessellation const& tess,
			bool cm_flag,
			SpatialDistribution const& density,
			SpatialDistribution const& pressure,
			SpatialDistribution const& xvelocity,
			SpatialDistribution const& yvelocity,
			EquationOfState const& eos) :
			tess_(tess), cm_flag_(cm_flag),
			density_(density),
			pressure_(pressure),
			xvelocity_(xvelocity),
			yvelocity_(yvelocity),
			eos_(eos) {}

		Primitive operator()(size_t n) const
		{
			return initialize_single_cell(tess_,
				static_cast<int>(n),
				cm_flag_,
				density_,
				pressure_,
				xvelocity_,
				yvelocity_,
				eos_);
		}

		size_t getLength(void) const
		{
			return static_cast<size_t>(tess_.GetPointNo());
		}

		~CellInitializer(void) {}

	private:
		Tessellation const& tess_;
		const bool cm_flag_;
		SpatialDistribution const& density_;
		SpatialDistribution const& pressure_;
		SpatialDistribution const& xvelocity_;
		SpatialDistribution const& yvelocity_;
		EquationOfState const& eos_;
	};
}

vector<Primitive> InitialiseCells
(SpatialDistribution const& density,
SpatialDistribution const& pressure,
SpatialDistribution const& xvelocity,
SpatialDistribution const& yvelocity,
EquationOfState const& eos,
Tessellation const& tess,
bool cm_value)
{
	return serial_generate
		(CellInitializer
		(tess, cm_value, density,
		pressure, xvelocity, yvelocity, eos));
}

namespace {
	class IntensiveInitializer : public Index2Member<Conserved>
	{
	public:

		IntensiveInitializer(vector<Primitive> const& cells) :
			cells_(cells) {}

		Conserved operator()(size_t n) const
		{
			return Primitive2Conserved(cells_[n]);
		}

		size_t getLength(void) const
		{
			return cells_.size();
		}

	private:
		vector<Primitive> const& cells_;
	};
}

vector<Conserved> CalcConservedIntensive
(vector<Primitive> const& cells)
{
	return serial_generate(IntensiveInitializer(cells));
}

namespace {

	class CellEdgesGetter : public Index2Member<Edge>
	{
	public:

		CellEdgesGetter(const Tessellation& tess, int n) :
			tess_(tess), edge_indices_(tess.GetCellEdges(n)) {}

		size_t getLength(void) const
		{
			return edge_indices_.size();
		}

		Edge operator()(size_t i) const
		{
			return tess_.GetEdge(edge_indices_[i]);
		}

	private:
		const Tessellation& tess_;
		const vector<int> edge_indices_;
	};

	class ExtensiveInitializer : public Index2Member<Conserved>
	{
	public:

		ExtensiveInitializer(const vector<Conserved>& intensive,
			const Tessellation& tess,
			const PhysicalGeometry& pg) :
			intensive_(intensive), tess_(tess), pg_(pg)  {}

		Conserved operator()(size_t n) const
		{
			return pg_.calcVolume(serial_generate(CellEdgesGetter(tess_, static_cast<int>(n))))*
				intensive_[n];
		}

		size_t getLength(void) const
		{
			return static_cast<size_t>(tess_.GetPointNo());
		}

	private:
		const vector<Conserved>& intensive_;
		const Tessellation& tess_;
		const PhysicalGeometry& pg_;
	};
}

vector<Conserved> CalcConservedExtensive
(const vector<Conserved>& cons_int,
const Tessellation& tess,
const PhysicalGeometry& pg)
{
	return serial_generate(ExtensiveInitializer(cons_int, tess, pg));
}

void CalcPointVelocities(Tessellation const& tessellation,
	vector<Primitive> const& cells,
	PointMotion& pointmotion,
	vector<Vector2D>& pointvelocity,
	double time,
	vector<CustomEvolution*> & cevolve,
	const vector<vector<double> >& tracers)
{
	pointvelocity = pointmotion.calcAllVelocities(tessellation, cells, time, cevolve, tracers);
}

double TimeStepForCell(Primitive const& cell, double width,
	vector<Vector2D> const& face_velocites)
{
	double max_fv = 0;
	for (size_t i = 0; i<face_velocites.size(); ++i)
		max_fv = max(max_fv, abs(face_velocites[i] - cell.Velocity));
	return width / (cell.SoundSpeed + max_fv);
}

namespace {
	double calc_boundary_face_velocity(HydroBoundaryConditions const& hbc,
		Edge const& edge,
		Tessellation const& tess,
		vector<Vector2D> const& face_velocities,
		Primitive const& cell,
		vector<Primitive> const& cells,
		double time,
		int i)
	{
		if (hbc.IsBoundary(edge, tess))
			return max(abs(face_velocities[static_cast<size_t>(i)] - cell.Velocity),
			abs(face_velocities[static_cast<size_t>(i)] -
			hbc.GetBoundaryPrimitive(edge,
			tess,
			cells,
			time).Velocity));
		else
			return abs(face_velocities[static_cast<size_t>(i)] - cell.Velocity);
	}
}

double TimeStepForCellBoundary
(Primitive const& cell,
vector<Primitive> const& cells,
double width,
vector<Vector2D> const& face_velocities,
Tessellation const& tess,
HydroBoundaryConditions const& hbc,
int index, double time)
{
	double max_fv = 0;
	const vector<int> edge_index = tess.GetCellEdges(index);
	for (int i = 0; i<static_cast<int>(face_velocities.size()); ++i)
	{
		const Edge& edge = tess.GetEdge(edge_index[static_cast<size_t>(i)]);
		const double temp = calc_boundary_face_velocity
			(hbc, edge, tess, face_velocities, cell,
			cells, time, i);
		max_fv = max(max_fv, temp);
	}
	return width / (cell.SoundSpeed + max_fv);
}

namespace {
	bool irrelevant_for_time_step(vector<CustomEvolution*> const& cevolve,
		int i)
	{
		if (!cevolve.empty())
			if (cevolve[static_cast<size_t>(i)] != 0)
				if (!cevolve[static_cast<size_t>(i)]->TimeStepRelevant())
					return true;
		return false;
	}
}

namespace {
	double calc_dt_temp(Tessellation const& tess,
		vector<Primitive> const& cells,
		vector<Vector2D> const& face_vel,
		HydroBoundaryConditions const& hbc,
		double time,
		int i)
	{
		if (!tess.NearBoundary(i))
			return TimeStepForCell(cells[static_cast<size_t>(i)], tess.GetWidth(i), face_vel);
		else
			return TimeStepForCellBoundary(cells[static_cast<size_t>(i)],
			cells,
			tess.GetWidth(i),
			face_vel,
			tess,
			hbc,
			i,
			time);
	}

	class FaceVelocityInitializer : public Index2Member<Vector2D>
	{
	public:

		FaceVelocityInitializer(vector<int> const& face_index,
			vector<Vector2D> const& face_velocity) :
			face_index_(face_index),
			face_velocity_(face_velocity) {}

		Vector2D operator()(size_t n) const
		{
			return face_velocity_[static_cast<size_t>(face_index_[n])];
		}

		size_t getLength(void) const
		{
			return face_index_.size();
		}

	private:
		vector<int> const& face_index_;
		vector<Vector2D> const& face_velocity_;
	};
}

double CalcTimeStep(Tessellation const& tessellation,
	vector<Primitive> const& cells,
	vector<Vector2D> const& facevelocity,
	HydroBoundaryConditions const& hbc,
	double time,
	vector<CustomEvolution*> const& evolve)
{
	bool first_time = true;
	double dt = 0;
	for (int i = 0; i<tessellation.GetPointNo(); i++)
	{
		if (irrelevant_for_time_step(evolve, i))
			continue;
		const vector<Vector2D> face_vel = serial_generate
			(FaceVelocityInitializer(tessellation.GetCellEdges(i),
			facevelocity));
		const double dt_temp = calc_dt_temp
			(tessellation, cells, face_vel,
			hbc, time, i);
		if (first_time)
		{
			first_time = false;
			dt = dt_temp;
		}
		else
			dt = min(dt_temp, dt);
	}
	return dt;
}

namespace {
#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
	__attribute__((noreturn))
#endif
		void update_conserved_extensive_error
		(int edge_index, int cell_number)
	{
		UniversalError eo("Error in UpdateConservedExtensive: cell and edge are not mutual neighbors");
		eo.AddEntry("edge number", edge_index);
		eo.AddEntry("cell number", cell_number);
		throw eo;
	}
}

void UpdateConservedExtensive
(Tessellation const& tessellation,
vector<Conserved> const& fluxes,
double dt,
vector<Conserved>& conserved_extensive,
HydroBoundaryConditions const& boundaryconditions,
vector<double> const& lengthes)
{
	for (int i = 0; i<tessellation.GetPointNo(); ++i){
		const vector<int> cell_edge_indices = tessellation.GetCellEdges(i);
		if (!boundaryconditions.IsGhostCell(i, tessellation)){
			for (int j = 0; j<static_cast<int>(cell_edge_indices.size()); ++j){
				const int edge_index = cell_edge_indices[static_cast<size_t>(j)];
				const Edge& edge = tessellation.GetEdge(edge_index);
				const Conserved delta = dt*lengthes[static_cast<size_t>(edge_index)] * fluxes[static_cast<size_t>(edge_index)];
				if (i == edge.neighbors.first)
					conserved_extensive[static_cast<size_t>(i)] -= delta;
				else if (i == edge.neighbors.second)
					conserved_extensive[static_cast<size_t>(i)] += delta;
				else
					update_conserved_extensive_error(edge_index, i);
			}
		}
	}
}

namespace {
	class NewPointPosition : public Index2Member<Vector2D>
	{
	public:

		NewPointPosition(Tessellation const& tess,
			vector<Vector2D> const& point_velocity,
			double dt) :
			tess_(tess),
			point_velocity_(point_velocity),
			dt_(dt) {}

		Vector2D operator()(size_t n) const
		{
			return tess_.GetMeshPoint(static_cast<int>(n)) + dt_*point_velocity_[n];
		}

		size_t getLength(void) const
		{
			return static_cast<size_t>(tess_.GetPointNo());
		}

	private:
		Tessellation const& tess_;
		vector<Vector2D> const& point_velocity_;
		const double dt_;
	};
}

void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
	double dt, Tessellation& tessellation, vector<Vector2D> oldpoints)
{
	if (oldpoints.empty())
	{
		tessellation.Update(serial_generate
			(NewPointPosition(tessellation,
			pointvelocity,
			dt)));
	}
	else
	{
		for (size_t i = 0; i<oldpoints.size(); ++i)
			oldpoints[i] += pointvelocity[i] * dt;
		tessellation.Update(oldpoints);
	}
}

#ifdef RICH_MPI
void MoveMeshPoints(vector<Vector2D> const& pointvelocity,
	double dt, Tessellation& tessellation,

	Tessellation const& vproc,
	vector<Vector2D> oldpoints)
{
	if (oldpoints.empty())
		tessellation.Update(serial_generate
		(NewPointPosition(tessellation,
		pointvelocity,
		dt)), vproc);
	else
	{
		for (size_t i = 0; i<oldpoints.size(); ++i)
			oldpoints[i] += pointvelocity[i] * dt;
		tessellation.Update(oldpoints, vproc);
	}
}
#endif // RICH_MPI

namespace {
	class IntensiveCalculator : public Index2Member<Conserved>
	{
	public:

		IntensiveCalculator(const Tessellation& tess,
			const vector<Conserved>& extensive,
			const PhysicalGeometry& pg) :
			tess_(tess), extensive_(extensive), pg_(pg) {}

		size_t getLength(void) const
		{
			return extensive_.size();
		}

		Conserved operator()(size_t i) const
		{
			return extensive_[i] /
				pg_.calcVolume(serial_generate(CellEdgesGetter(tess_, static_cast<int>(i))));
		}

	private:
		const Tessellation& tess_;
		const vector<Conserved>& extensive_;
		const PhysicalGeometry& pg_;
	};
}

vector<Conserved> calc_conserved_intensive
(const Tessellation& tess,
const vector<Conserved>& extensive,
const PhysicalGeometry& pg)
{
	return serial_generate(IntensiveCalculator(tess, extensive, pg));
}

void UpdateConservedIntensive(Tessellation const& tessellation,
	vector<Conserved> const& conservedextensive,
	vector<Conserved>& conservedintensive)
{
	conservedintensive.resize(conservedextensive.size());
	for (int i = 0; i<tessellation.GetPointNo(); i++){
		conservedintensive[static_cast<size_t>(i)] = conservedextensive[static_cast<size_t>(i)] /
			tessellation.GetVolume(i);
	}
}

namespace {

	std::pair<Conserved, bool> calc_safe_conserved(Conserved const& raw,
		bool density_floor,
		double min_density,
		double min_pressure,
		Primitive const& old,
		EquationOfState const& eos)
	{
		std::pair<Conserved, bool> res;
		res.first = raw;
		if (density_floor)
		{
			if (res.first.Mass < min_density)
			{
				res.first.Mass = min_density;
				res.first.Momentum = old.Velocity*min_density;
				const double kinetic_energy = 0.5*pow(abs(res.first.Momentum / res.first.Mass), 2);
				res.first.Energy = res.first.Mass*kinetic_energy
					+ res.first.Mass*eos.dp2e(res.first.Mass, min_pressure);
				res.second = true;
			}
			const double kinetic_energy = 0.5*pow(abs(res.first.Momentum / res.first.Mass), 2);
			const double thermal_energy = res.first.Energy / res.first.Mass - kinetic_energy;
			const double pressure = eos.de2p(res.first.Mass, thermal_energy);
			if (pressure < min_pressure || res.second)
			{
				res.first.Energy = res.first.Mass*kinetic_energy
					+ res.first.Mass*eos.dp2e(res.first.Mass, min_pressure);
				res.second = true;
			}
		}
		return res;
	}

#if defined(__clang__) || defined(__GNUC__) || defined(__GNUG__)
	__attribute__((noreturn))
#endif
		void update_primitives_rethrow(int cell_index,
		UniversalError& eo)
	{
		eo.AddEntry("UpdatePrimitive data starts here", 0);
		eo.AddEntry("cell index", static_cast<double>(cell_index));
		throw eo;
	}

	std::pair<Primitive, bool> regular_cell_evolve(Conserved const& intensive,
		bool density_floor,
		double min_density,
		double min_pressure,
		Primitive const& old,
		EquationOfState const& eos)
	{
		const std::pair<Conserved, bool> temp = calc_safe_conserved
			(intensive, density_floor, min_density,
			min_pressure, old, eos);
		return std::pair<Primitive, bool>(Conserved2Primitive(temp.first, eos), temp.second);
	}
}

vector<bool> UpdatePrimitives
(vector<Conserved> const& conservedintensive,
EquationOfState const& eos, vector<Primitive>& cells,
vector<CustomEvolution*> const& CellsEvolve, vector<Primitive> &old_cells,
bool densityfloor, double densitymin, double pressuremin, Tessellation const&
tess, double time, vector<vector<double> > const& tracers)
{
	cells.resize(conservedintensive.size());
	vector<bool> bres(cells.size(), false);
	for (int i = 0; i < tess.GetPointNo(); i++)
	{
		try
		{
			if (CellsEvolve[static_cast<size_t>(i)] == 0 || CellsEvolve[static_cast<size_t>(i)]->DensityFloorRelevant())
			{
				Primitive old_cell = densityfloor ? old_cells[static_cast<size_t>(i)] : cells[static_cast<size_t>(i)];
				std::pair<Primitive, bool > res = regular_cell_evolve
					(conservedintensive[static_cast<size_t>(i)], densityfloor,
					densitymin, pressuremin, old_cell, eos);
				if (CellsEvolve[static_cast<size_t>(i)] != 0 && !res.second)
					cells[static_cast<size_t>(i)] = CellsEvolve[static_cast<size_t>(i)]->UpdatePrimitive
					(conservedintensive, eos, old_cells, i, tess, time, tracers);
				else
					cells[static_cast<size_t>(i)] = res.first;
				bres[static_cast<size_t>(i)] = res.second;
			}
			else
				cells[static_cast<size_t>(i)] = CellsEvolve[static_cast<size_t>(i)]->UpdatePrimitive
				(conservedintensive,
				eos, old_cells, i, tess, time, tracers);
		}
		catch (UniversalError& eo)
		{
			eo.AddEntry("x momentum per unit volume", conservedintensive[static_cast<size_t>(i)].Momentum.x);
			eo.AddEntry("y momentum per unit volume", conservedintensive[static_cast<size_t>(i)].Momentum.y);
			eo.AddEntry("thermal energy per unit mass", conservedintensive[static_cast<size_t>(i)].Energy);
			eo.AddEntry("Cell volume", tess.GetVolume(i));
			eo.AddEntry("Cell x location", tess.GetMeshPoint(i).x);
			eo.AddEntry("Cell y location", tess.GetMeshPoint(i).y);
#ifdef RICH_MPI
			eo.AddEntry("Error in CPU", static_cast<double>(get_mpi_rank()));
#endif
			update_primitives_rethrow(i, eo);
		}
	}
	return bres;
}

Primitive RotatePrimitive(Vector2D const& normaldir,
	Vector2D const& paraldir,
	Primitive const& p)
{
	Primitive res = p;
	res.Velocity.Set(Projection(p.Velocity, normaldir),
		Projection(p.Velocity, paraldir));
	return res;
}

Conserved RotateFluxBack(Conserved const& c,
	Vector2D const& normaldir,
	Vector2D const& paraldir)
{
	Conserved res = c;
	res.Momentum = c.Momentum.x*normaldir / abs(normaldir) +
		c.Momentum.y*paraldir / abs(paraldir);
	return res;
}

Conserved FluxInBulk(Vector2D const& normaldir,
	Vector2D const& paraldir,
	Primitive const& left,
	Primitive const& right,
	Vector2D const& edge_velocity,
	RiemannSolver const& rs)
{
	const Primitive rotated_left = RotatePrimitive(normaldir, paraldir, left);
	const Primitive rotated_right = RotatePrimitive(normaldir, paraldir, right);
	const double normal_speed = Projection(edge_velocity, normaldir);
	const Conserved res = rs.Solve(rotated_left, rotated_right, normal_speed);
	return RotateFluxBack(res, normaldir, paraldir);
}

namespace {
	Conserved calc_single_flux_in_bulk(Tessellation const& tess,
		Edge const& edge,
		SpatialReconstruction const& interpolation,
		vector<Primitive> const& cells,
		Vector2D const& face_velocity,
		RiemannSolver const& rs,
		double dt)
	{
		const Vector2D normal_dir =
			tess.GetMeshPoint(edge.neighbors.second) -
			tess.GetMeshPoint(edge.neighbors.first);

		const Vector2D paral_dir =
			edge.vertices.second - edge.vertices.first;

		const Primitive left = interpolation.Interpolate
			(tess, cells, dt, edge, 0, InBulk, face_velocity);

		const Primitive right = interpolation.Interpolate
			(tess, cells, dt, edge, 1, InBulk, face_velocity);

		return FluxInBulk(normal_dir, paral_dir,
			left, right,
			face_velocity, rs);
	}

	int choose_special_cell_index
		(vector<CustomEvolution*> const& ce_list,
		const CustomEvolutionManager& cem,
		int n0, int n1)
	{
		if (ce_list[static_cast<size_t>(n0)] && ce_list[static_cast<size_t>(n1)]){
			if (ce_list[static_cast<size_t>(n0)] == ce_list[static_cast<size_t>(n1)])
				return n0;
			const size_t priority_0 = cem.getIndex(ce_list[static_cast<size_t>(n0)]);
			const size_t priority_1 = cem.getIndex(ce_list[static_cast<size_t>(n1)]);
			assert(priority_0 != priority_1 &&
				"Two methods have the same priority");
			return priority_0 > priority_1 ?
			n0 : n1;
		}
		else if (!ce_list[static_cast<size_t>(n0)] && !ce_list[static_cast<size_t>(n1)])
			throw UniversalError("Error in choose_special_cell_index: Both sides are regular cells");
		else if (ce_list[static_cast<size_t>(n0)] && !ce_list[static_cast<size_t>(n1)])
			return n0;
		else if (!ce_list[static_cast<size_t>(n0)] && ce_list[static_cast<size_t>(n1)])
			return n1;
		else
			throw UniversalError("Error in choose_special_cell_index: Something has gone terribly wrong if you've gotten here");
	}

	/*
	void calc_fluxes_rethrow(UniversalError& eo,
	int edge_index,
	Tessellation const& tess)
	{
	eo.AddEntry("Error in CalcFlux",0);
	eo.AddEntry("edge index",edge_index);
	const Edge& edge = tess.GetEdge(edge_index);
	eo.AddEntry("edge x1 location",edge.vertices.first.x);
	eo.AddEntry("edge y1 location",edge.vertices.first.y);
	eo.AddEntry("edge x2 location",edge.vertices.second.x);
	eo.AddEntry("edge y2 location",edge.vertices.second.y);
	throw eo;
	}
	*/
}

namespace {
	class InterpolationRelevancy : public Index2Member<bool>
	{
	public:

		InterpolationRelevancy(const vector<CustomEvolution*>& ce) :
			ce_(ce) {}

		size_t getLength(void) const
		{
			return ce_.size();
		}

		bool operator()(size_t i) const
		{
			if (ce_[i])
				return ce_[i]->isRelevantToInterpolation();
			else
				return true;
		}

	private:
		const vector<CustomEvolution*>& ce_;
	};

	class FluxCalculator : public Index2Member<Conserved>
	{
	public:

		FluxCalculator(SpatialReconstruction& interp,
			const Tessellation& tess,
			const vector<Primitive>& cells,
			const vector<vector<double> >& tracers,
			const vector<CustomEvolution*>& ce,
			const CustomEvolutionManager& cem,
			double time, double dt,
			const HydroBoundaryConditions& hbc,
			const vector<Vector2D>& fv,
			const RiemannSolver& rs) :
			interp_(interp),
			tess_(tess),
			cells_(cells),
			tracers_(tracers),
			hbc_(hbc),
			ce_(ce),
			cem_(cem),
			fv_(fv),
			rs_(rs),
			time_(time),
			dt_(dt)
		{
			interp_.Prepare(tess, cells, tracers,
				serial_generate(InterpolationRelevancy(ce)),
				dt, time);
#ifndef RICH_MPI
			PeriodicGradExchange(interp_.GetGradients(),
				tess_.GetDuplicatedPoints(), tess_.GetTotalPointNumber());
#else
			SendRecvGrad(interp_.GetGradients(), tess_.GetDuplicatedPoints(),
				tess_.GetDuplicatedProcs(), tess_.GetGhostIndeces(),
				tess_.GetTotalPointNumber());
#endif
		}

		size_t getLength(void) const
		{
			return static_cast<size_t>(tess_.GetTotalSidesNumber());
		}

		Conserved operator()(size_t i) const
		{
			const Edge& edge = tess_.GetEdge(static_cast<int>(i));
			const int n0 = edge.neighbors.first;
			const int n1 = edge.neighbors.second;
			if (!hbc_.IsBoundary(edge, tess_)){
				if (!ce_[static_cast<size_t>(n0)] && !ce_[static_cast<size_t>(n1)])
					return calc_single_flux_in_bulk(tess_, edge, interp_, cells_, fv_[i], rs_, dt_);
				else{
					const int ns = choose_special_cell_index(ce_, cem_, n0, n1);
					return ce_[static_cast<size_t>(ns)]->CalcFlux(tess_, cells_, dt_, interp_, edge, fv_[i], rs_, ns, hbc_, time_, tracers_);
				}
			}
			else
				return hbc_.CalcFlux(tess_, cells_, fv_[i], edge, interp_, dt_, time_);
		}

	private:
		SpatialReconstruction& interp_;
		const Tessellation& tess_;
		const vector<Primitive>& cells_;
		const vector<vector<double> >& tracers_;
		const HydroBoundaryConditions& hbc_;
		const vector<CustomEvolution*>& ce_;
		const CustomEvolutionManager& cem_;
		const vector<Vector2D>& fv_;
		const RiemannSolver& rs_;
		const double time_;
		const double dt_;
	};
}

vector<Conserved> calc_fluxes
(Tessellation const& tessellation,
vector<Primitive> const& cells,
double dt,
double time,
SpatialReconstruction& interpolation,
vector<Vector2D> const& facevelocity,
HydroBoundaryConditions const& boundaryconditions,
RiemannSolver const& rs,
vector<CustomEvolution*> const& CellsEvolve,
CustomEvolutionManager const& cem,
vector<vector<double> > const& tracers)
{
	return serial_generate
		(FluxCalculator(interpolation,
		tessellation,
		cells,
		tracers,
		CellsEvolve,
		cem,
		time, dt,
		boundaryconditions,
		facevelocity,
		rs));
}

void ExternalForceContribution
(Tessellation const& tess,
const PhysicalGeometry& pg,
vector<Primitive> const& cells,
SourceTerm& force,
double t,
double dt,
vector<Conserved>& conserved_extensive,
HydroBoundaryConditions const& hbc, vector<Conserved> const& fluxes,
vector<Vector2D> const& point_velocity, vector<double> &g, bool coldflows_flag,
vector<vector<double> > &tracers_extensive, vector<double> const& lengthes)
{
	vector<double> dtracer;
	if (!tracers_extensive.empty())
		dtracer.assign(tracers_extensive[0].size(), 0);
	if (coldflows_flag)
		g.resize(static_cast<size_t>(tess.GetPointNo()));
	for (int i = 0; i<tess.GetPointNo(); ++i)
	{
		if (!hbc.IsGhostCell(i, tess))
		{
			const Conserved cons(force.Calculate
				(tess, pg, cells, i,
				fluxes, point_velocity, hbc, tracers_extensive,
				dtracer, lengthes, t, dt));
			conserved_extensive[static_cast<size_t>(i)] += dt*cons;
			if (!tracers_extensive.empty())
			{
				tracers_extensive[static_cast<size_t>(i)] = tracers_extensive[static_cast<size_t>(i)] + dt*dtracer;
			}
			if (coldflows_flag)
				g[static_cast<size_t>(i)] = abs(cons.Momentum) / (cells[static_cast<size_t>(i)].Density*tess.GetVolume(i));
		}
	}
}

vector<Vector2D> calc_point_velocities
(Tessellation const& tess,
vector<Primitive> const& cells,
PointMotion& point_motion,
double time, vector<CustomEvolution*> & cevolve,
const vector<vector<double> >& tracers)
{
	return point_motion.calcAllVelocities(tess, cells, time, cevolve, tracers);
}

vector<Vector2D> get_all_mesh_points
(Tessellation const& tess)
{
	vector<Vector2D> res(static_cast<size_t>(tess.GetPointNo()));
	for (int i = 0; i<static_cast<int>(tess.GetPointNo()); ++i)
		res[static_cast<size_t>(i)] = tess.GetMeshPoint(i);
	return res;
}

vector<CustomEvolution*> convert_indices_to_custom_evolution
(const CustomEvolutionManager& cem, const vector<size_t>& indices)
{
	vector<CustomEvolution*> res(indices.size());
	for (size_t i = 0; i<res.size(); ++i)
		res[i] = cem.getFunction(indices[i]);
	return res;
}

namespace
{
	void FirstOrderLargeChanges(vector<Conserved> const& old_extensive, vector<Conserved> &extensive,
		vector<Conserved> const& dextensive, vector<vector<double> > const&old_tracers,
		vector<vector<double> > &tracers, vector<vector<double> > const& dtracer)
	{
		const double factor = 2;
		for (size_t i = 0; i < extensive.size(); ++i)
		{
			const double e_ratio = std::abs(dextensive[i].Energy) / old_extensive[i].Energy;
			if (std::abs(dextensive[i].Mass) > factor*old_extensive[i].Mass || e_ratio>factor)
			{
				extensive[i] = old_extensive[i] + dextensive[i];
				if (!tracers.empty())
				{
					for (size_t j = 0; j < tracers[i].size(); ++j)
						tracers[i][j] = old_tracers[i][j] + dtracer[i][j];
				}
			}
		}
	}
}

double TimeAdvance2mid(Tessellation& tess,
#ifdef RICH_MPI
Tessellation& proctess,
#endif
vector<Primitive> &cells,PointMotion& point_motion, HydroBoundaryConditions const& hbc,
SpatialReconstruction& interpolation, RiemannSolver const& rs,
EquationOfState const& eos, SourceTerm& force, double time, double cfl,
double endtime, vector<vector<double> > &tracers,double dt_external,
vector<size_t>& custom_evolution_indices,const CustomEvolutionManager& custom_evolution_manager,
const PhysicalGeometry& pg,
#ifdef RICH_MPI
ProcessorUpdate *procupdate,
#endif
bool traceflag, bool coldflows_flag,double as, double bs, bool densityfloor, double densitymin,
double pressuremin, bool EntropyCalc)
{
	vector<Primitive> old_cells = cells;
	vector<vector<double> > old_trace_intensive = tracers;
	// create ghost points if needed
#ifndef RICH_MPI
	PeriodicUpdateCells(cells, tracers, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());

	PeriodicVelocityExchange(tess.GetAllCM(), tess.GetDuplicatedPoints(), tess.GetTotalPointNumber());
#else
	SendRecvHydro(cells, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), eos, tess.GetGhostIndeces(), tess.GetTotalPointNumber());

	SendRecvVelocity(tess.GetAllCM(), tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif
	vector<CustomEvolution*> CellsEvolve = convert_indices_to_custom_evolution(
		custom_evolution_manager, custom_evolution_indices);
	vector<Vector2D> oldpoints = tess.GetMeshPoints();
	oldpoints.resize(static_cast<size_t>(tess.GetPointNo()));

	//do half time step
	vector<Vector2D> point_velocities = calc_point_velocities
		(tess, cells, point_motion, time, CellsEvolve, tracers);

#ifndef RICH_MPI
	PeriodicVelocityExchange(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());
#else
	SendRecvVelocity(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif

	vector<Vector2D> edge_velocities = tess.calc_edge_velocities
		(&hbc, point_velocities, time);

	double dt = determine_time_step
		(cfl*CalcTimeStep(tess, cells, edge_velocities, hbc, time, CellsEvolve),
		dt_external, time, endtime);

	point_motion.ApplyFix(tess, cells, time, CellsEvolve, tracers, 0.5*dt, point_velocities);
#ifndef RICH_MPI
	PeriodicVelocityExchange(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());
#else
	SendRecvVelocity(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif
	edge_velocities = tess.calc_edge_velocities(&hbc, point_velocities, time);
	dt = determine_time_step(cfl*CalcTimeStep(tess, cells, edge_velocities, hbc, time, CellsEvolve),
		dt_external, time, endtime);

	// Entropy and tracers evolution
	vector<double> g, Ek, Ef;
	vector<char> shockedcells;
	if (coldflows_flag)
	{
		const int n = tess.GetPointNo();
		for (int i = 0; i<n; ++i)
			tracers[static_cast<size_t>(i)][0] = eos.dp2s(cells[static_cast<size_t>(i)].Density, cells[static_cast<size_t>(i)].Pressure);
#ifdef RICH_MPI
		SendRecvTracers(tracers, tess.GetDuplicatedPoints(),
			tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
		if (tracers.size() != cells.size())
		{
			UniversalError eo("Tracers and cells have different length at first half time step");
			eo.AddEntry("CPU rank", get_mpi_rank());
			throw eo;
		}
#endif
		shockedcells.resize(static_cast<size_t>(n));
		for (int i = 0; i<n; ++i)
			shockedcells[static_cast<size_t>(i)] = IsShockedCell(tess, i, cells, hbc, time) ? 1 : 0;
	}
#ifdef RICH_MPI
	if (traceflag&&!coldflows_flag)
	{
		SendRecvTracers(tracers, tess.GetDuplicatedPoints(),
			tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
		if (tracers.size() != cells.size())
		{
			UniversalError eo("Tracers and cells have different length at first half time step");
			eo.AddEntry("CPU rank", get_mpi_rank());
			throw eo;
		}
	}
#endif

	vector<Conserved> fluxes = calc_fluxes
		(tess, cells, 0.5*dt, time, interpolation,
		edge_velocities, hbc, rs, CellsEvolve,
		custom_evolution_manager, tracers);

	vector<double> lengths = serial_generate(EdgeLengthCalculator(tess, pg));

	vector<vector<double> > old_trace;
	vector<vector<double> > tracer_extensive;
	if (traceflag)
	{
		vector<vector<double> > trace_change;
		trace_change = CalcTraceChange
			(tracers, cells, tess, fluxes, 0.5*dt, hbc,
			interpolation, time, CellsEvolve,
			custom_evolution_manager,
			edge_velocities, lengths);
		MakeTracerExtensive(tracers, tess, cells, tracer_extensive);
		old_trace = tracer_extensive;
		UpdateTracerExtensive(tracer_extensive, trace_change, CellsEvolve, cells,
			tess, time);
	}

	vector<Conserved> intensive = CalcConservedIntensive(cells);

	vector<Conserved> extensive = CalcConservedExtensive
		(intensive, tess, pg);

	// Save extensive variables of beginning of time step
	vector<Conserved> old_extensive = extensive;

	UpdateConservedExtensive(tess, fluxes, 0.5*dt,
		extensive, hbc, lengths);

	ExternalForceContribution
		(tess, pg, cells, force, time, 0.5*dt,
		extensive, hbc, fluxes, point_velocities, g, coldflows_flag, tracers, lengths);

	vector<Conserved> dextensive = 2 * (extensive - old_extensive);
	vector<vector<double> > dtracer = 2 * (tracer_extensive - old_trace);

	if (coldflows_flag)
	{
		Ek = GetMaxKineticEnergy(tess, cells, CellsEvolve);
		Ef = GetForceEnergy(tess, g);
	}

#ifndef RICH_MPI
	MoveMeshPoints(point_velocities, 0.5*dt, tess);
#else
	MoveMeshPoints(point_velocities, 0.5*dt, tess, proctess);
#endif

#ifdef RICH_MPI
	vector<Conserved> ptoadd;
	vector<vector<double> > ttoadd;
	vector<size_t> ctemp(custom_evolution_indices);
	vector<size_t> ctoadd;
	SendRecvExtensive(extensive, tracer_extensive, custom_evolution_indices,
		tess.GetSentPoints(), tess.GetSentProcs(), ptoadd, ttoadd,
		ctoadd);

	vector<Conserved> dctoadd;
	vector<vector<double> > dttoadd;
	vector<size_t> ctoadd2;
	SendRecvExtensive(dextensive, dtracer, custom_evolution_indices,
		tess.GetSentPoints(), tess.GetSentProcs(), dctoadd, dttoadd,
		ctoadd2);

	vector<size_t> ctemp2(custom_evolution_indices);
	KeepLocalPoints(dextensive, dtracer, ctemp2,
		tess.GetSelfPoint());
	if (!dctoadd.empty())
		dextensive.insert(dextensive.end(), dctoadd.begin(), dctoadd.end());

	if (!dttoadd.empty())
		dtracer.insert(dtracer.end(), dttoadd.begin(), dttoadd.end());

	KeepLocalPoints(extensive, tracer_extensive, custom_evolution_indices,
		tess.GetSelfPoint());

	if (!ptoadd.empty())
	{
		extensive.insert(extensive.end(), ptoadd.begin(),
			ptoadd.end());
	}

	if (!ttoadd.empty())
	{
		tracer_extensive.insert(tracer_extensive.end(), ttoadd.begin(), ttoadd.end());
	}

	if (!ctoadd.empty())
	{
		custom_evolution_indices.insert(custom_evolution_indices.end(), ctoadd.begin(),
			ctoadd.end());
	}

	SendRecvExtensive(old_extensive, old_trace, ctemp,
		tess.GetSentPoints(), tess.GetSentProcs(), ptoadd, ttoadd,
		ctoadd);

	KeepLocalPoints(old_extensive, old_trace, ctemp, tess.GetSelfPoint());

	if (!ptoadd.empty())
	{
		old_extensive.insert(old_extensive.end(), ptoadd.begin(),
			ptoadd.end());
	}

	if (!ttoadd.empty())
	{
		old_trace.insert(old_trace.end(), ttoadd.begin(), ttoadd.end());
	}

	vector<Vector2D> vtoadd;
	SendRecvOldVector2D(oldpoints, tess.GetSentPoints(),
		tess.GetSentProcs(), vtoadd);
	oldpoints = VectorValues(oldpoints, tess.GetSelfPoint());
	if (!vtoadd.empty())
		oldpoints.insert(oldpoints.end(), vtoadd.begin(), vtoadd.end());

	vtoadd.clear();
	SendRecvOldVector2D(point_velocities, tess.GetSentPoints(),
		tess.GetSentProcs(), vtoadd);
	point_velocities = VectorValues(point_velocities, tess.GetSelfPoint());
	if (!vtoadd.empty())
		point_velocities.insert(point_velocities.end(), point_velocities.begin(), point_velocities.end());

	CellsEvolve =
		convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	if (coldflows_flag)
	{
		vector<char> btoadd;
		SendRecvShockedStatus(shockedcells, tess.GetSentPoints(), tess.GetSentProcs(),
			btoadd);
		shockedcells = VectorValues(shockedcells, tess.GetSelfPoint());
		if (!btoadd.empty())
			shockedcells.insert(shockedcells.end(), btoadd.begin(), btoadd.end());
		vector<double> Ekadd;
		SendRecvVectorDouble(Ek, tess.GetSentPoints(), tess.GetSentProcs(),
			Ekadd);
		Ek = VectorValues(Ek, tess.GetSelfPoint());
		if (!Ekadd.empty())
			Ek.insert(Ek.end(), Ekadd.begin(), Ekadd.end());
		SendRecvVectorDouble(Ef, tess.GetSentPoints(), tess.GetSentProcs(),
			Ekadd);
		Ef = VectorValues(Ef, tess.GetSelfPoint());
		if (!Ekadd.empty())
			Ef.insert(Ef.end(), Ekadd.begin(), Ekadd.end());
	}
	if (densityfloor)
	{
		vector<Primitive> pmtoadd;
		ttoadd.clear();
		SendRecvPrimitive(old_cells, old_trace_intensive, tess.GetSentPoints(), tess.GetSentProcs(),
			pmtoadd, ttoadd, eos);
		KeepLocalPoints(old_cells, old_trace_intensive, tess.GetSelfPoint());
		if (!pmtoadd.empty())
			old_cells.insert(old_cells.end(), pmtoadd.begin(), pmtoadd.end());
		if (!ttoadd.empty())
			old_trace_intensive.insert(old_trace_intensive.begin(), ttoadd.begin(), ttoadd.end());
	}
#endif

	UpdateConservedIntensive(tess, extensive, intensive);

	if (coldflows_flag)
		FixPressure(intensive, tracer_extensive, eos, Ek, Ef, as, bs, CellsEvolve,
		tess,/*extensive,*/shockedcells, densityfloor);

	cells.resize(static_cast<size_t>(tess.GetPointNo()));
	vector<bool> min_density_on = UpdatePrimitives(intensive, eos, cells, CellsEvolve, old_cells, densityfloor,
		densitymin, pressuremin, tess, time + 0.5*dt, tracers);
	if (traceflag)
		MakeTracerIntensive(tracers, tracer_extensive, tess, cells, pg, min_density_on, old_trace, CellsEvolve);

	// End half step

	// create ghost points if needed
#ifndef RICH_MPI
	PeriodicUpdateCells(cells, tracers, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());

	PeriodicVelocityExchange(tess.GetAllCM(), tess.GetDuplicatedPoints(), tess.GetTotalPointNumber());
#else
	SendRecvHydro(cells, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), eos, tess.GetGhostIndeces(), tess.GetTotalPointNumber());

	SendRecvVelocity(tess.GetAllCM(), tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif

	CellsEvolve = convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

#ifndef RICH_MPI
	PeriodicVelocityExchange(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());
#else
	SendRecvVelocity(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif

	edge_velocities = tess.calc_edge_velocities(&hbc, point_velocities, time + 0.5*dt);

	if (coldflows_flag&&EntropyCalc)
	{
		const int n = tess.GetPointNo();
		for (size_t i = 0; i<static_cast<size_t>(n); ++i)
			tracers[i][0] = eos.dp2s(cells[i].Density, cells[i].Pressure);
#ifdef RICH_MPI
		SendRecvTracers(tracers, tess.GetDuplicatedPoints(),
			tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
		if (tracers.size() != cells.size())
		{
			UniversalError eo("Tracers and cells have different length at second half time step");
			eo.AddEntry("CPU rank", get_mpi_rank());
			throw eo;
		}

#endif
	}

#ifdef RICH_MPI
	if (traceflag && (!coldflows_flag || !EntropyCalc))
	{
		SendRecvTracers(tracers, tess.GetDuplicatedPoints(),
			tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
		if (tracers.size() != cells.size())
		{
			UniversalError eo("Tracers and cells have different length at second half time step");
			eo.AddEntry("CPU rank", get_mpi_rank());
			throw eo;
		}
	}
#endif
	fluxes = calc_fluxes
		(tess, cells, dt, time, interpolation,
		edge_velocities, hbc, rs, CellsEvolve,
		custom_evolution_manager, tracers);

	if (coldflows_flag)
	{
		const int n = tess.GetPointNo();
		shockedcells.resize(static_cast<size_t>(n));
		for (int i = 0; i<n; ++i)
			shockedcells[static_cast<size_t>(i)] = IsShockedCell(tess, i, cells, hbc, time) ? 1 : 0;
	}

	lengths = serial_generate(EdgeLengthCalculator(tess, pg));

	extensive = old_extensive;

	if (traceflag)
	{
		old_trace_intensive = tracers;
		vector<vector<double> > trace_change;
		trace_change = CalcTraceChange
			(tracers, cells, tess, fluxes, dt, hbc,
			interpolation, time, CellsEvolve,
			custom_evolution_manager,
			edge_velocities, lengths);
		tracer_extensive = old_trace;
		UpdateTracerExtensive(old_trace, trace_change, CellsEvolve, cells,
			tess, time);
	}

	UpdateConservedExtensive(tess, fluxes, dt,
		old_extensive, hbc, lengths);

	ExternalForceContribution
		(tess, pg, cells, force, time + 0.5*dt, dt,
		old_extensive, hbc, fluxes, point_velocities, g, coldflows_flag, tracers, lengths);

	if (coldflows_flag)
	{
		Ek = GetMaxKineticEnergy(tess, cells, CellsEvolve);
		Ef = GetForceEnergy(tess, g);
	}

	FirstOrderLargeChanges(extensive, old_extensive, dextensive, tracer_extensive, old_trace, dtracer);

#ifndef RICH_MPI
	MoveMeshPoints(point_velocities, dt, tess, oldpoints);
#else
	if (procupdate != 0)
		procupdate->Update(proctess, tess);
	MoveMeshPoints(point_velocities, dt, tess, proctess, oldpoints);
#endif

#ifdef RICH_MPI

	SendRecvExtensive(old_extensive, old_trace, custom_evolution_indices,
		tess.GetSentPoints(), tess.GetSentProcs(), ptoadd, ttoadd,
		ctoadd);

	KeepLocalPoints(old_extensive, old_trace, custom_evolution_indices,
		tess.GetSelfPoint());

	if (!ptoadd.empty())
	{
		old_extensive.insert(old_extensive.end(), ptoadd.begin(),
			ptoadd.end());
	}

	if (!ttoadd.empty())
	{
		old_trace.insert(old_trace.end(), ttoadd.begin(), ttoadd.end());
	}

	if (!ctoadd.empty())
	{
		custom_evolution_indices.insert(custom_evolution_indices.end(), ctoadd.begin(),
			ctoadd.end());
	}

	CellsEvolve =
		convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	if (coldflows_flag)
	{
		vector<char> btoadd;
		SendRecvShockedStatus(shockedcells, tess.GetSentPoints(), tess.GetSentProcs(),
			btoadd);
		shockedcells = VectorValues(shockedcells, tess.GetSelfPoint());
		if (!btoadd.empty())
			shockedcells.insert(shockedcells.end(), btoadd.begin(), btoadd.end());
		vector<double> Ekadd;
		SendRecvVectorDouble(Ek, tess.GetSentPoints(), tess.GetSentProcs(),
			Ekadd);
		Ek = VectorValues(Ek, tess.GetSelfPoint());
		if (!Ekadd.empty())
			Ek.insert(Ek.end(), Ekadd.begin(), Ekadd.end());
		SendRecvVectorDouble(Ef, tess.GetSentPoints(), tess.GetSentProcs(),
			Ekadd);
		Ef = VectorValues(Ef, tess.GetSelfPoint());
		if (!Ekadd.empty())
			Ef.insert(Ef.end(), Ekadd.begin(), Ekadd.end());
	}

	if (densityfloor)
	{
		old_cells = cells;
		vector<Primitive> pmtoadd;
		ttoadd.clear();
		SendRecvPrimitive(old_cells, old_trace_intensive, tess.GetSentPoints(), tess.GetSentProcs(),
			pmtoadd, ttoadd, eos);
		KeepLocalPoints(old_cells, old_trace_intensive, tess.GetSelfPoint());
		if (!pmtoadd.empty())
			old_cells.insert(old_cells.end(), pmtoadd.begin(), pmtoadd.end());
		if (!ttoadd.empty())
			old_trace_intensive.insert(old_trace_intensive.begin(), ttoadd.begin(), ttoadd.end());
	}
#endif

	UpdateConservedIntensive(tess, old_extensive, intensive);
	//intensive = calc_conserved_intensive(tess, old_extensive, pg);

	if (coldflows_flag)
		FixPressure(intensive, old_trace, eos, Ek, Ef, as, bs, CellsEvolve, tess,
		/*extensive,*/shockedcells, densityfloor);

	min_density_on = UpdatePrimitives
		(intensive, eos, cells, CellsEvolve, old_cells, densityfloor,
		densitymin, pressuremin, tess, time + dt, tracers);

	if (traceflag)
	{
		MakeTracerIntensive(tracers, old_trace, tess, cells, pg, min_density_on, old_trace_intensive, CellsEvolve);
	}
	return dt;
}

vector<Primitive> make_eos_consistent
(vector<Primitive> const& vp,
EquationOfState const& eos)
{
	vector<Primitive> res = vp;
	for (int i = 0; i<static_cast<int>(vp.size()); ++i)
		res[static_cast<size_t>(i)] = make_eos_consistent(vp[static_cast<size_t>(i)], eos);
	return res;
}

namespace {
	/*
	class ConditionalPlusMinus: public BinaryOperation<double>
	{
	public:

	ConditionalPlusMinus(bool flag):
	flag_(flag) {}

	double operator()(double const& x, double const& y) const
	{
	if(flag_)
	return x+y;
	else
	return x-y;
	}

	private:
	const bool flag_;
	};
	*/

	class ScalarMultiply : public UnaryOperation<double>
	{
	public:

		ScalarMultiply(double scalar) :
			scalar_(scalar) {}

		double operator()(double const& x) const
		{
			return x*scalar_;
		}

	private:
		const double scalar_;
	};
}

namespace {
	class TracerFluxCalculator : public Index2Member<vector<double> >
	{
	public:

		TracerFluxCalculator(const Tessellation& tess,
			const vector<Conserved> & fluxes,
			const vector<CustomEvolution*>& cev,
			const CustomEvolutionManager& cem,
			const vector<vector<double> >& tracers,
			const vector<Primitive>& cells,
			const HydroBoundaryConditions& hbc,
			const double dt,
			const double time,
			const SpatialReconstruction& interp,
			const vector<Vector2D>& edge_velocities,
			const vector<double>& lengths) :
			tess_(tess), fluxes_(fluxes), cev_(cev), cem_(cem),
			tracers_(tracers), cells_(cells), hbc_(hbc), dt_(dt),
			time_(time), interp_(interp),
			edge_velocities_(edge_velocities), lengths_(lengths) {}

		size_t getLength(void) const
		{
			return static_cast<size_t>(tess_.GetTotalSidesNumber());
		}

		vector<double> operator()(size_t i) const
		{
			const Edge& edge = tess_.GetEdge(static_cast<int>(i));
			const double dm = fluxes_[i].Mass;
			const int n1 = edge.neighbors.second;
			const int n0 = edge.neighbors.first;
			if (hbc_.IsBoundary(edge, tess_)){
				const bool b0 = hbc_.IsGhostCell(n0, tess_);
				const bool b1 = hbc_.IsGhostCell(n1, tess_);
				assert(b0 != b1 && "One cell must be a normal cell, and the other must be a ghost cell");
				const int rci = b1 ? n0 : n1;
				return hbc_.CalcTracerFlux
					(tess_, cells_, tracers_, dm, edge, rci, dt_, time_, interp_,
					edge_velocities_[i]);
			}
			else{
				if (cev_[static_cast<size_t>(n0)] || cev_[static_cast<size_t>(n1)])
					return cev_[static_cast<size_t>(choose_special_cell_index(cev_, cem_, n0, n1))]->CalcTracerFlux
					(tess_, cells_, tracers_, dm, edge, static_cast<int>(i), dt_,
					time_, interp_, edge_velocities_[i]);
				else
					return apply_to_each_term
					(interp_.interpolateTracers
					(tess_, cells_, tracers_, dt_, edge, dm<0, InBulk,
					edge_velocities_[i]),
					ScalarMultiply(dt_*dm*lengths_[i]));
			}
		}

	private:
		const Tessellation& tess_;
		const vector<Conserved>& fluxes_;
		const vector<CustomEvolution*>& cev_;
		const CustomEvolutionManager& cem_;
		const vector<vector<double> >& tracers_;
		const vector<Primitive>& cells_;
		const HydroBoundaryConditions& hbc_;
		const double dt_;
		const double time_;
		const SpatialReconstruction& interp_;
		const vector<Vector2D>& edge_velocities_;
		const vector<double>& lengths_;
	};
}

void really_update_extensive_tracers
(vector<vector<double> >& extensive_tracers,
const vector<vector<double> >& tracers,
const vector<Primitive>& cells,
const Tessellation& tess,
const vector<Conserved>& fluxes,
double time, double dt,
const HydroBoundaryConditions& hbc,
const SpatialReconstruction& interp,
const vector<CustomEvolution*>& ce,
const CustomEvolutionManager& cem,
const vector<Vector2D>& fv,
const vector<double>& lengths)
{
	const vector<vector<double> >& tracer_change =
		CalcTraceChange(tracers, cells, tess, fluxes, dt, hbc,
		interp, time, ce, cem, fv, lengths);
	UpdateTracerExtensive
		(extensive_tracers, tracer_change,
		ce, cells, tess, time);
}

vector<vector<double> > CalcTraceChange
(vector<vector<double> > const& old_trace,
vector<Primitive> const& cells,
Tessellation const& tess, vector<Conserved> const& fluxes, double dt,
HydroBoundaryConditions const& hbc,
SpatialReconstruction const& interp,
double time, vector<CustomEvolution*> const& CellsEvolve,
CustomEvolutionManager const& cem,
vector<Vector2D> const& edge_velocities,
vector<double> const& lengths)
{
	if (old_trace.empty())
		return vector<vector<double> >();

	const vector<vector<double> >& tracer_fluxes =
		serial_generate(TracerFluxCalculator
		(tess, fluxes, CellsEvolve, cem,
		old_trace, cells, hbc, dt, time,
		interp, edge_velocities, lengths));
	vector<vector<double> > res(old_trace.size(),
		vector<double>(old_trace[0].size(), 0));
	for (int i = 0; i<tess.GetTotalSidesNumber(); ++i){
		for (size_t j = 0; j<res[0].size(); ++j){
			if (!hbc.IsGhostCell(tess.GetEdge(i).neighbors.first, tess))
				res[static_cast<size_t>(tess.GetEdge(i).neighbors.first)][j] -= tracer_fluxes[static_cast<size_t>(i)][j];
			if (!hbc.IsGhostCell(tess.GetEdge(i).neighbors.second, tess))
				res[static_cast<size_t>(tess.GetEdge(i).neighbors.second)][j] += tracer_fluxes[static_cast<size_t>(i)][j];
		}
	}
	return res;
}

vector<double> GetMaxKineticEnergy(Tessellation const& tess, vector<Primitive> const&
	cells, vector<CustomEvolution*> const& /*customevolve*/)
{
	const int n = tess.GetPointNo();
	vector<double> res;
	res.resize(static_cast<size_t>(n));
	for (int j = 0; j<n; ++j)
	{
		vector<int> neightemp = tess.GetNeighbors(j);
		vector<int> neigh;
		for (size_t i = 0; i<neightemp.size(); ++i)
			if (neightemp[i] >= 0)
				neigh.push_back(neightemp[i]);
		double e = pow(abs(cells[static_cast<size_t>(j)].Velocity - cells[static_cast<size_t>(neigh[0])].Velocity), 2);
		for (int i = 1; i<static_cast<int>(neigh.size()); ++i)
		{// This could be made much faster by writing the expression implicitly
			e = max(e, pow(abs(cells[static_cast<size_t>(j)].Velocity - cells[static_cast<size_t>(neigh[static_cast<size_t>(i)])].Velocity), 2));
		}
		res[static_cast<size_t>(j)] = 0.5*e;
	}
	return res;
}

vector<double> GetForceEnergy(Tessellation const& tess,
	vector<double> const& g)
{
	vector<double> res;
	int n = int(g.size());
	res.resize(static_cast<size_t>(n));
	for (int i = 0; i<n; ++i)
		res[static_cast<size_t>(i)] = g[static_cast<size_t>(i)] * tess.GetWidth(i);
	return res;
}

void FixPressure(vector<Conserved> &intensive, vector<vector<double> > const& entropy,
	EquationOfState const& eos, vector<double> const& Ek,
	vector<double> const& Ef, double as, double bs, vector<CustomEvolution*>
	const& customevolve, Tessellation const& tess,//vector<Conserved> &extensive,
	vector<char> const& shockedcells, bool densityfloor)
{
	int n = tess.GetPointNo();
	double Et, Ek2;
	double temp;
	for (int i = 0; i<n; ++i)
	{
		if (customevolve[static_cast<size_t>(i)] == 0 || customevolve[static_cast<size_t>(i)]->TimeStepRelevant())
		{
			if (intensive[static_cast<size_t>(i)].Mass < 0)
				continue;
			//Make intensive
			temp = entropy[static_cast<size_t>(i)][0] / (tess.GetVolume(i)*intensive[static_cast<size_t>(i)].Mass);
			Ek2 = 0.5*pow(abs(intensive[static_cast<size_t>(i)].Momentum) / intensive[static_cast<size_t>(i)].Mass, 2);
			Et = intensive[static_cast<size_t>(i)].Energy / intensive[static_cast<size_t>(i)].Mass - Ek2;
			if ((Et<as*Ek[static_cast<size_t>(i)]) || (Et<bs*Ef[static_cast<size_t>(i)]))
			{
				if ((shockedcells[static_cast<size_t>(i)] == 0) || Et<0)
				{
					Et = eos.dp2e(intensive[static_cast<size_t>(i)].Mass,
						eos.sd2p(temp, intensive[static_cast<size_t>(i)].Mass));
					if (Et<0 && !densityfloor)
					{
						UniversalError eo("Negative thermal enegry");
						eo.AddEntry("Cell index", i);
						eo.AddEntry("Thermal energy", Et);
						eo.AddEntry("ShockedStatus", shockedcells[static_cast<size_t>(i)]);
						eo.AddEntry("Extensive entropy", entropy[static_cast<size_t>(i)][0]);
						eo.AddEntry("The density", intensive[static_cast<size_t>(i)].Mass);
						eo.AddEntry("X cor", tess.GetMeshPoint(i).x);
						eo.AddEntry("Y cor", tess.GetMeshPoint(i).y);
						throw eo;
					}
					intensive[static_cast<size_t>(i)].Energy = intensive[static_cast<size_t>(i)].Mass*(Et + Ek2);
					//extensive[i].Energy=tess.GetVolume(i)*intensive[i].Energy;
				}
			}
		}
	}
}

bool NearBoundary(int index, Tessellation const& tess,
	vector<CustomEvolution*> const& /*customevolve*/)
{
	vector<int> neigh = tess.GetNeighbors(index);
	int n = int(neigh.size());
	for (int i = 0; i<n; ++i)
	{
		if (neigh[static_cast<size_t>(i)]<0)
			return true;
		/*if(customevolve[neigh[i]]!=0)
		return true;*/
	}
	return false;
}

namespace {
	vector<double> scalar_mult(const vector<double>& v,
		double s)
	{
		if (v.empty())
			return vector<double>();
		vector<double> res(v.size());
		for (size_t i = 0; i<v.size(); ++i)
			res[i] = s*v[i];
		return res;
	}
}

namespace {
	class ExtensiveTracerCalculator : public Index2Member<vector<double> >
	{
	public:

		ExtensiveTracerCalculator(const vector<vector<double> >& tracers,
			const Tessellation& tess,
			const vector<Primitive>& cells,
			const PhysicalGeometry& pg) :
			tracers_(tracers), tess_(tess), cells_(cells), pg_(pg) {}

		size_t getLength(void) const
		{
			if (tracers_.empty())
				return 0;
			else
				return static_cast<size_t>(tess_.GetPointNo());
		}

		vector<double> operator()(size_t i) const
		{
			return scalar_mult
				(tracers_[i],
				pg_.calcVolume(serial_generate(CellEdgesGetter(tess_, static_cast<int>(i))))*
				cells_[i].Density);
		}

	private:
		const vector<vector<double> >& tracers_;
		const Tessellation& tess_;
		const vector<Primitive>& cells_;
		const PhysicalGeometry& pg_;
	};
}

vector<vector<double> > calc_extensive_tracer
(const vector<vector<double> >& intensive_tracer,
const Tessellation& tess,
const vector<Primitive>& cells,
const PhysicalGeometry& pg)
{
	return serial_generate(ExtensiveTracerCalculator(intensive_tracer,
		tess,
		cells,
		pg));
}

void MakeTracerExtensive(vector<vector<double> > const &tracer,
	Tessellation const& tess,
	vector<Primitive> const& cells,
	vector<vector<double> > &result)
{
	const size_t n = static_cast<size_t>(tess.GetPointNo());
	result.resize(n);
	for (size_t i = 0; i<n; ++i)
		result[i] = scalar_mult(tracer[i],
		tess.GetVolume(static_cast<int>(i))*cells[i].Density);
}

namespace {
	vector<double> scalar_div(const vector<double>& v,
		const double s)
	{
		vector<double> res(v.size());
		for (size_t i = 0; i<res.size(); ++i)
			res[i] = v[i] / s;
		return res;
	}

	class IntensiveTracerCalculator : public Index2Member<vector<double> >
	{
	public:

		IntensiveTracerCalculator(const vector<vector<double> >& extensive, const Tessellation& tess,
			const vector<Primitive>& cells, const PhysicalGeometry& pg,
			const vector<vector<double> > &old_trace_intensive, const vector<bool> &min_density,
			const vector<CustomEvolution*> &cevolve) :
			extensive_(extensive), tess_(tess), cells_(cells), pg_(pg),
			old_trace_intensive_(old_trace_intensive), min_density_(min_density), cevolve_(cevolve){}

		size_t getLength(void) const
		{
			if (extensive_.empty())
				return 0;
			else
				return static_cast<size_t>(tess_.GetPointNo());
		}

		vector<double> operator()(size_t i) const
		{
			if (min_density_[i] && !cevolve_[i])
				return old_trace_intensive_[i];
			else
			{
				const double mass = cells_[i].Density*
					pg_.calcVolume(serial_generate(CellEdgesGetter(tess_, static_cast<int>(i))));
				return scalar_div(extensive_[i], mass);
			}
		}

	private:
		const vector<vector<double> >& extensive_;
		const Tessellation& tess_;
		const vector<Primitive>& cells_;
		const PhysicalGeometry& pg_;
		const vector<vector<double> > &old_trace_intensive_;
		const vector<bool> &min_density_;
		const vector<CustomEvolution*> &cevolve_;
	};
}

void MakeTracerIntensive(vector<vector<double> > &tracer,
	const vector<vector<double> >& extensive,
	const Tessellation& tess,
	const vector<Primitive>& cells,
	const PhysicalGeometry& pg, vector<bool> const& min_density_on, vector<vector<double> > const& old_trace,
	vector<CustomEvolution*> const& cevolve)
{
	tracer = serial_generate(IntensiveTracerCalculator(extensive, tess, cells, pg, old_trace, min_density_on, cevolve));
}

void UpdateTracerExtensive(vector<vector<double> > &tracerextensive,
	vector<vector<double> > const& tracerchange, vector<CustomEvolution*> const&
	CellsEvolve, vector<Primitive> const& cells, Tessellation const& tess,
	double time)
{
	for (size_t i = 0; i<tracerextensive.size(); ++i)
		if (CellsEvolve[i])
			tracerextensive[i] = CellsEvolve[i]->UpdateTracer
			(static_cast<int>(i), tracerextensive, tracerchange, cells, tess, time);
		else
			for (size_t j = 0; j<tracerextensive[i].size(); ++j)
				tracerextensive[i][j] += tracerchange[i][j];
}

Vector2D TracerResetCalc
(double alpha, SpatialDistribution const& originalD,
SpatialDistribution const& originalP, SpatialDistribution const& originalVx,
SpatialDistribution const& originalVy, vector<SpatialDistribution const*> const& originalTracers, vector<Primitive> &cells,
Tessellation const& tess, vector<vector<double> > &tracer,
int tracerindex, EquationOfState const& eos, vector<CustomEvolution*>
const& cevolve, bool coldflows)
{
	Vector2D res(0, 0);
	const int n = tess.GetPointNo();
	if (n<1)
		return res;
	Vector2D velocity;
	if (tracer.empty())
		return res;
	if (tracerindex >= static_cast<int>(tracer[0].size()) || tracerindex<0)
		throw UniversalError("Error in tracerReset, wrong dimension for tracer");
	for (int i = 0; i<n; ++i)
	{
		bool customforce = false;
		if (cevolve[static_cast<size_t>(i)])
			customforce = cevolve[static_cast<size_t>(i)]->ShouldForceTracerReset();
		if ((tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)]<alpha) || customforce)
		{
			velocity.Set(originalVx(tess.GetCellCM(i)),
				originalVy(tess.GetCellCM(i)));
			double mass = cells[static_cast<size_t>(i)].Density*tess.GetVolume(i);
			double tmass = mass*tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)];
			cells[static_cast<size_t>(i)] = CalcPrimitive(originalD(tess.GetCellCM(i)),
				originalP(tess.GetCellCM(i)), velocity, eos);
			if (tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)]<0)
				tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)] = 0;
			double toadd = cells[static_cast<size_t>(i)].Density*tess.GetVolume(i) - mass;
			res.x += toadd;
			for (size_t j = 0; j<tracer[static_cast<size_t>(i)].size(); ++j)
				if (coldflows&&j == 0)
					tracer[static_cast<size_t>(i)][j] = eos.dp2s(cells[static_cast<size_t>(i)].Density, cells[static_cast<size_t>(i)].Pressure);
				else
					if (static_cast<int>(j) != tracerindex)
						tracer[static_cast<size_t>(i)][j] = originalTracers[j]->operator()(tess.GetCellCM(i));
			res.y += tmass - cells[static_cast<size_t>(i)].Density*tess.GetVolume(i)*tracer[static_cast<size_t>(i)][static_cast<size_t>(tracerindex)];
		}
	}
	return res;
}

void GetPointToRemove(Tessellation const& tess, Vector2D const& point,
	double R, vector<int> & PointToRemove, int Inner)
{
	int n = tess.GetPointNo();
	PointToRemove.clear();
	for (int i = Inner; i<n; ++i)
	{
		// Check if point is completly engulfed
		bool test = true;
		vector<int> neigh = tess.GetNeighbors(i);
		for (int j = 0; j<static_cast<int>(neigh.size()); ++j)
			if (neigh[static_cast<size_t>(j)] >= Inner)
				test = false;
		// Is point inside a radius?
		if (abs(point - tess.GetMeshPoint(i))<R || test)
			PointToRemove.push_back(i);
	}
	return;
}

namespace {
	Vector2D GetReflectedPoint(Tessellation const& tess, int point,
		Edge const& edge)
	{
		Vector2D MeshPoint = tess.GetMeshPoint(point);
		Vector2D par = edge.vertices.second - edge.vertices.first;
		if (fabs(par.x)>fabs(par.y)){
			// We are in the x direction
			if (MeshPoint.y>edge.vertices.first.y)
				MeshPoint.y -= MeshPoint.y - edge.vertices.first.y;
			else
				MeshPoint.y += edge.vertices.first.y - MeshPoint.y;
		}
		else{
			// We are in the y direction
			if (MeshPoint.x>edge.vertices.first.x)
				MeshPoint.x -= MeshPoint.x - edge.vertices.first.x;
			else
				MeshPoint.x += edge.vertices.first.x - MeshPoint.x;
		}
		return MeshPoint;
	}
}

bool IsShockedCell(Tessellation const& tess, int index,
	vector<Primitive> const& cells,
	HydroBoundaryConditions const& hbc,
	double time)
{
	double DivV = 0;
	vector<int> edges_loc = tess.GetCellEdges(index);
	vector<Edge> edges(edges_loc.size());
	for (size_t i = 0; i<edges.size(); ++i)
		edges[i] = tess.GetEdge(edges_loc[i]);

	// Calculate gradiant by gauss theorem
	double vx_i = cells[static_cast<size_t>(index)].Velocity.x;
	double vy_i = cells[static_cast<size_t>(index)].Velocity.y;
	Vector2D center = tess.GetMeshPoint(index);
	Vector2D c_ij, r_ij;
	for (size_t i = 0; i<edges.size(); i++)
	{
		if (hbc.IsBoundary(edges[i], tess))
		{
			Primitive other2 = hbc.GetBoundaryPrimitive(edges[i], tess, cells,
				time);
			if ((edges[i].neighbors.first == -1) || (edges[i].neighbors.second == -1))
			{
				c_ij = CalcCentroid(edges[i]) - 0.5*(GetReflectedPoint
					(tess, index, edges[i]) + center);
				r_ij = center - GetReflectedPoint(tess, index, edges[i]);
			}
			else
			{
				const int other_point = get_other_index(edges[i], index);
				c_ij = CalcCentroid(edges[i]) - 0.5*(
					tess.GetMeshPoint(other_point) + center);
				r_ij = center - tess.GetMeshPoint(other_point);
			}
			double rij_1 = 1 / abs(r_ij);
			double vx_j = other2.Velocity.x;
			double vy_j = other2.Velocity.y;
			DivV += edges[i].GetLength()*((vx_j - vx_i)*c_ij.x - 0.5*
				(vx_i + vx_j)*r_ij.x + (vy_j - vy_i)*c_ij.y - 0.5*
				(vy_i + vy_j)*r_ij.y)*rij_1;
		}
		else
		{
			const int other = get_other_index(edges[i], index);
			c_ij = CalcCentroid(edges[i]) -
				0.5*(tess.GetMeshPoint(other) + center);
			r_ij = center - tess.GetMeshPoint(other);
			double rij_1 = 1 / abs(r_ij);
			double vx_j = cells[static_cast<size_t>(other)].Velocity.x;
			double vy_j = cells[static_cast<size_t>(other)].Velocity.y;
			DivV += edges[i].GetLength()*((vx_j - vx_i)*c_ij.x - 0.5*
				(vx_i + vx_j)*r_ij.x + (vy_j - vy_i)*c_ij.y - 0.5*
				(vy_i + vy_j)*r_ij.y)*rij_1;
		}
	}
	if (DivV<-0.1*tess.GetWidth(index)*cells[static_cast<size_t>(index)].SoundSpeed)
		return true;
	else
		return false;
}

void FixAdvection(vector<Conserved>& extensive,
	vector<Conserved> const& intensive, Tessellation const& tessold,
	Tessellation const& tessnew, vector<Vector2D> const& facevelocity,
	double dt, vector<Vector2D> const& /*pointvelocity*/)
{
	int n = tessold.GetTotalSidesNumber();
	int npoints = tessold.GetPointNo();
	vector<double> Rold(static_cast<size_t>(npoints)), Rnew(static_cast<size_t>(npoints));
	vector<vector<Vector2D> > pold(static_cast<size_t>(npoints)), pnew(static_cast<size_t>(npoints));
	for (int i = 0; i<npoints; ++i)
	{
		Rold[static_cast<size_t>(i)] = tessold.GetWidth(i);
		Rnew[static_cast<size_t>(i)] = tessnew.GetWidth(i);
		ConvexHull(pold[static_cast<size_t>(i)], &tessold, i);
		ConvexHull(pnew[static_cast<size_t>(i)], &tessnew, i);
	}

	PolygonOverlap polyoverlap;
	double eps = 1e-7;
	for (int i = 0; i<n; ++i)
	{
		Edge const& edge = tessold.GetEdge(i);
		int n0 = edge.neighbors.first;
		int n1 = edge.neighbors.second;
		if (n0<0 || n1<0)
			continue;
		Vector2D norm(tessold.GetMeshPoint(n1) - tessold.GetMeshPoint(n0));
		norm = norm / abs(norm);
		norm = norm*edge.GetLength();
		double dv_dt = ScalarProd(facevelocity[static_cast<size_t>(i)], norm)*dt;
		/*		vector<Vector2D> poly0,poly1;
		ConvexHull(poly0,&tessold,tessold.GetOriginalIndex(n0));
		ConvexHull(poly1,&tessnew,tessold.GetOriginalIndex(n1));
		if(n0>=npoints)
		{
		const Vector2D diff(tessold.GetMeshPoint(tessold.GetOriginalIndex(n0))
		-tessold.GetMeshPoint(n0));
		int N=static_cast<int>(poly0.size());
		for(int j=0;j<N;++j)
		poly0[j]-=diff;
		}
		if(n1>=npoints)
		{
		const Vector2D diff(tessnew.GetMeshPoint(tessnew.GetOriginalIndex(n1))
		-tessnew.GetMeshPoint(n1));
		int N=static_cast<int>(poly1.size());
		for(int j=0;j<N;++j)
		poly1[j]-=diff;
		}*/
		double real_dv1 = polyoverlap.polygon_overlap_area
			(pold[static_cast<size_t>(tessold.GetOriginalIndex(n0))],
			pnew[static_cast<size_t>(tessold.GetOriginalIndex(n1))],
			Rold[static_cast<size_t>(tessold.GetOriginalIndex(n0))] * eps,
			Rnew[static_cast<size_t>(tessold.GetOriginalIndex(n1))] * eps);
		double real_dv0 = polyoverlap.polygon_overlap_area
			(pnew[static_cast<size_t>(tessold.GetOriginalIndex(n0))],
			pold[static_cast<size_t>(tessold.GetOriginalIndex(n1))],
			Rnew[static_cast<size_t>(tessold.GetOriginalIndex(n0))] * eps,
			Rold[static_cast<size_t>(tessold.GetOriginalIndex(n1))] * eps);

		if (dv_dt>0)
		{
			if (n0<npoints)
			{
				extensive[static_cast<size_t>(n0)] += (real_dv0 - dv_dt)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
				extensive[static_cast<size_t>(n0)] -= real_dv1*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
			}
			if (n1<npoints)
			{
				extensive[static_cast<size_t>(n1)] += (dv_dt - real_dv0)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
				extensive[static_cast<size_t>(n1)] += real_dv1*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
			}
		}
		else
		{
			if (n0<npoints)
			{
				extensive[static_cast<size_t>(n0)] -= (real_dv1 + dv_dt)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
				extensive[static_cast<size_t>(n0)] += real_dv0*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
			}
			if (n1<npoints)
			{
				extensive[static_cast<size_t>(n1)] -= (-dv_dt - real_dv1)*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n0))];
				extensive[static_cast<size_t>(n1)] -= real_dv0*intensive[static_cast<size_t>(tessold.GetOriginalIndex(n1))];
			}
		}
	}
}

double determine_time_step(double hydro_time_step,
	double external_dt,
	double current_time,
	double end_time)
{
	double dt = hydro_time_step;
	if (external_dt>0)
		dt = std::min(external_dt, dt);
	if (end_time>0)
		dt = std::min(end_time - current_time, dt);

#ifdef RICH_MPI
	double dt_temp = dt;
	MPI_Reduce(&dt_temp, &dt, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

	return dt;
}

double TimeAdvance2midClip(Tessellation& tess, vector<Primitive> &cells,
	PointMotion& point_motion, HydroBoundaryConditions const& hbc,
	SpatialReconstruction& interpolation, RiemannSolver const& rs,
	EquationOfState const& eos, SourceTerm& force, double time, double cfl,
	double endtime, vector<vector<double> > &tracers, double dt_external,
	vector<size_t>& custom_evolution_indices,
	const CustomEvolutionManager& custom_evolution_manager,
	const PhysicalGeometry& pg, bool traceflag, bool coldflows_flag, double as,
	double bs, bool densityfloor, double densitymin, double pressuremin,
	bool EntropyCalc, OuterBoundary const& outer)
{
	vector<Primitive> old_cells = cells;
	vector<vector<double> > old_trace_intensive = tracers;
	// create ghost points if needed
	PeriodicUpdateCells(cells, tracers, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());

	PeriodicVelocityExchange(tess.GetAllCM(), tess.GetDuplicatedPoints(), tess.GetTotalPointNumber());

	vector<CustomEvolution*> CellsEvolve = convert_indices_to_custom_evolution(
		custom_evolution_manager, custom_evolution_indices);

	vector<Vector2D> oldpoints = tess.GetMeshPoints();
	oldpoints.resize(static_cast<size_t>(tess.GetPointNo()));

	//do half time step
	vector<Vector2D> point_velocities = calc_point_velocities
		(tess, cells, point_motion, time, CellsEvolve, tracers);

	PeriodicVelocityExchange(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());

	vector<Vector2D> edge_velocities = tess.calc_edge_velocities
		(&hbc, point_velocities, time);

	double dt = determine_time_step
		(cfl*CalcTimeStep(tess, cells, edge_velocities, hbc, time, CellsEvolve),
		dt_external, time, endtime);

	point_motion.ApplyFix(tess, cells, time, CellsEvolve, tracers, 0.5*dt, point_velocities);
	PeriodicVelocityExchange(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());
	edge_velocities = tess.calc_edge_velocities(&hbc, point_velocities, time);
	dt = determine_time_step(cfl*CalcTimeStep(tess, cells, edge_velocities, hbc, time, CellsEvolve),
		dt_external, time, endtime);

	// Entropy and tracers evolution
	vector<double> g, Ek, Ef;
	vector<char> shockedcells;
	if (coldflows_flag)
	{
		const int n = tess.GetPointNo();
		for (int i = 0; i<n; ++i)
		  tracers[static_cast<size_t>(i)][0] = eos.dp2s(cells[static_cast<size_t>(i)].Density, cells[static_cast<size_t>(i)].Pressure);
		shockedcells.resize(static_cast<size_t>(n));
		for (int i = 0; i<n; ++i)
		  shockedcells[static_cast<size_t>(i)] = IsShockedCell(tess, i, cells, hbc, time) ? 1 : 0;
	}

	vector<Conserved> fluxes = calc_fluxes
		(tess, cells, 0.5*dt, time, interpolation,
		edge_velocities, hbc, rs, CellsEvolve,
		custom_evolution_manager, tracers);

	vector<double> lengths = serial_generate(EdgeLengthCalculator(tess, pg));


	vector<vector<double> > old_trace;
	vector<vector<double> > tracer_extensive;
	if (traceflag)
	{
		vector<vector<double> > trace_change;
		trace_change = CalcTraceChange(tracers, cells, tess, fluxes, 0.5*dt, hbc,
			interpolation, time, CellsEvolve, custom_evolution_manager,
			edge_velocities, lengths);
		MakeTracerExtensive(tracers, tess, cells, tracer_extensive);
		old_trace = tracer_extensive;
		UpdateTracerExtensive(tracer_extensive, trace_change, CellsEvolve, cells,
			tess, time);
	}

	vector<Conserved> intensive = CalcConservedIntensive(cells);

	vector<Conserved> extensive = CalcConservedExtensive
		(intensive, tess, pg);

	// Save extensive variables of beginning of time step
	vector<Conserved> old_extensive = extensive;

	boost::scoped_ptr<Tessellation> oldtess(tess.clone());

	MoveMeshPoints(point_velocities, 0.5*dt, tess);

	vector<vector<double>  > ttoadd;
	vector<Conserved> toadd = FluxFix2(*oldtess, *oldtess, tess, point_velocities, 0.5*dt, cells, fluxes, edge_velocities,
		outer, eos, hbc, tracers, ttoadd);
	extensive = extensive + toadd;
	if (!tracers.empty())
		tracer_extensive = tracer_extensive + ttoadd;

	UpdateConservedExtensive(*oldtess, fluxes, 0.5*dt,
		extensive, hbc, lengths);

	ExternalForceContribution(*oldtess, pg, cells, force, time, 0.5*dt,
		extensive, hbc, fluxes, point_velocities, g, coldflows_flag, tracers, lengths);

	vector<Conserved> dextensive = 2 * (extensive - old_extensive);
	vector<vector<double> > dtracer = 2 * (tracer_extensive - old_trace);


	if (coldflows_flag)
	{
		Ek = GetMaxKineticEnergy(tess, cells, CellsEvolve);
		Ef = GetForceEnergy(tess, g);
	}

	UpdateConservedIntensive(tess, extensive, intensive);

	if (coldflows_flag)
		FixPressure(intensive, tracer_extensive, eos, Ek, Ef, as, bs, CellsEvolve,
		tess,/*extensive,*/shockedcells, densityfloor);

	cells.resize(static_cast<size_t>(tess.GetPointNo()));
	vector<bool> min_density_on = UpdatePrimitives(intensive, eos, cells, CellsEvolve, cells, densityfloor,
		densitymin, pressuremin, tess, time + 0.5*dt, tracers);
	if (traceflag)
	{
		MakeTracerIntensive(tracers, tracer_extensive, tess, cells, pg, min_density_on, old_trace, CellsEvolve);
	}

	// End half step

	// create ghost points if needed
	PeriodicUpdateCells(cells, tracers, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());

	PeriodicVelocityExchange(tess.GetAllCM(), tess.GetDuplicatedPoints(), tess.GetTotalPointNumber());

	CellsEvolve = convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	PeriodicVelocityExchange(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());

	edge_velocities = tess.calc_edge_velocities(&hbc, point_velocities, time);

	if (coldflows_flag&&EntropyCalc)
	{
		const int n = tess.GetPointNo();
		for (size_t i = 0; i<static_cast<size_t>(n); ++i)
			tracers[i][0] = eos.dp2s(cells[i].Density, cells[i].Pressure);
	}

	fluxes = calc_fluxes
		(tess, cells, dt, time, interpolation,
		edge_velocities, hbc, rs, CellsEvolve,
		custom_evolution_manager, tracers);

	if (coldflows_flag)
	{
		const int n = tess.GetPointNo();
		shockedcells.resize(static_cast<size_t>(n));
		for (int i = 0; i<n; ++i)
		  shockedcells[static_cast<size_t>(i)] = IsShockedCell(tess, i, cells, hbc, time) ? 1 : 0;
	}

	lengths = serial_generate(EdgeLengthCalculator(tess, pg));

	extensive = old_extensive;
	tracer_extensive = old_trace;

	if (traceflag)
	{
		vector<vector<double> > trace_change;
		trace_change = CalcTraceChange
			(tracers, cells, tess, fluxes, dt, hbc,
			interpolation, time, CellsEvolve,
			custom_evolution_manager,
			edge_velocities, lengths);
		UpdateTracerExtensive(old_trace, trace_change, CellsEvolve, cells,
			tess, time);
	}

	boost::scoped_ptr<Tessellation> midtess(tess.clone());

	MoveMeshPoints(point_velocities, dt, tess, oldpoints);

	old_extensive = old_extensive + FluxFix2(*oldtess, *midtess, tess, point_velocities, dt, cells, fluxes, edge_velocities,
		outer, eos, hbc, tracers, ttoadd);
	if (!tracers.empty())
		old_trace = old_trace + ttoadd;

	UpdateConservedExtensive(*midtess, fluxes, dt,
		old_extensive, hbc, lengths);

	ExternalForceContribution
		(*midtess, pg, cells, force, time + 0.5*dt, dt,
		old_extensive, hbc, fluxes, point_velocities, g, coldflows_flag, tracers, lengths);

	if (coldflows_flag)
	{
		Ek = GetMaxKineticEnergy(tess, cells, CellsEvolve);
		Ef = GetForceEnergy(tess, g);
	}

	FirstOrderLargeChanges(extensive, old_extensive, dextensive, tracer_extensive, old_trace, dtracer);


	UpdateConservedIntensive(tess, old_extensive, intensive);

	if (coldflows_flag)
		FixPressure(intensive, old_trace, eos, Ek, Ef, as, bs, CellsEvolve, tess,
		/*extensive,*/shockedcells, densityfloor);

	min_density_on = UpdatePrimitives
		(intensive, eos, cells, CellsEvolve, cells, densityfloor,
		densitymin, pressuremin, tess, time + dt, tracers);

	if (traceflag)
	{
		MakeTracerIntensive(tracers, old_trace, tess, cells, pg, min_density_on, old_trace_intensive, CellsEvolve);
	}
	return dt;
}

double TimeAdvance2New(Tessellation& tess,
#ifdef RICH_MPI
	Tessellation& proctess,
#endif
	vector<Primitive> &cells, PointMotion& point_motion, HydroBoundaryConditions const& hbc, SpatialReconstruction& interpolation, RiemannSolver const& rs,
	EquationOfState const& eos, SourceTerm& force, double time, double cfl, double endtime, vector<vector<double> > &tracers,
	double dt_external, vector<size_t>& custom_evolution_indices, const CustomEvolutionManager& custom_evolution_manager,
	const PhysicalGeometry& pg,
#ifdef RICH_MPI
	ProcessorUpdate *procupdate,
#endif
	bool traceflag, bool coldflows_flag, double as, double bs, bool densityfloor, double densitymin, double pressuremin, bool EntropyCalc)
{
	// create ghost points if needed
#ifndef RICH_MPI
	PeriodicUpdateCells(cells, tracers, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());

	PeriodicVelocityExchange(tess.GetAllCM(), tess.GetDuplicatedPoints(), tess.GetTotalPointNumber());
#else
	SendRecvHydro(cells, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), eos, tess.GetGhostIndeces(), tess.GetTotalPointNumber());

	SendRecvVelocity(tess.GetAllCM(), tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif
	vector<CustomEvolution*> CellsEvolve = convert_indices_to_custom_evolution(
		custom_evolution_manager, custom_evolution_indices);

	//do first time step
	vector<Vector2D> point_velocities = calc_point_velocities(tess, cells, point_motion, time, CellsEvolve, tracers);

#ifndef RICH_MPI
	PeriodicVelocityExchange(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());
#else
	SendRecvVelocity(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif

	vector<Vector2D> edge_velocities = tess.calc_edge_velocities(&hbc, point_velocities, time);

	double dt = determine_time_step(cfl*CalcTimeStep(tess, cells, edge_velocities, hbc, time, CellsEvolve),
		dt_external, time, endtime);

	point_motion.ApplyFix(tess, cells, time, CellsEvolve, tracers, dt, point_velocities);
#ifndef RICH_MPI
	PeriodicVelocityExchange(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());
#else
	SendRecvVelocity(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif
	edge_velocities = tess.calc_edge_velocities(&hbc, point_velocities, time);

	// Entropy and tracers evolution
	vector<double> g, Ek, Ef;
	vector<char> shockedcells;
	if (coldflows_flag)
	{
		const int n = tess.GetPointNo();
		for (int i = 0; i<n; ++i)
			tracers[static_cast<size_t>(i)][0] = eos.dp2s(cells[static_cast<size_t>(i)].Density, cells[static_cast<size_t>(i)].Pressure);
#ifdef RICH_MPI
		SendRecvTracers(tracers, tess.GetDuplicatedPoints(),
			tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
		if (tracers.size() != cells.size())
		{
			UniversalError eo("Tracers and cells have different length at first half time step");
			eo.AddEntry("CPU rank", get_mpi_rank());
			throw eo;
		}
#endif
		shockedcells.resize(static_cast<size_t>(n));
		for (int i = 0; i<n; ++i)
			shockedcells[static_cast<size_t>(i)] = IsShockedCell(tess, i, cells, hbc, time) ? 1 : 0;
	}
#ifdef RICH_MPI
	if (traceflag&&!coldflows_flag)
	{
		SendRecvTracers(tracers, tess.GetDuplicatedPoints(),
			tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
		if (tracers.size() != cells.size())
		{
			UniversalError eo("Tracers and cells have different length at first half time step");
			eo.AddEntry("CPU rank", get_mpi_rank());
			throw eo;
		}
	}
#endif

	vector<Conserved> fluxes = calc_fluxes(tess, cells, dt, time, interpolation, edge_velocities, hbc, rs, CellsEvolve,
		custom_evolution_manager, tracers);
	vector<double> lengths = serial_generate(EdgeLengthCalculator(tess, pg));
	vector<vector<double> > old_trace_extensive;
	vector<vector<double> > tracer_extensive;
	if (traceflag)
	{
		vector<vector<double> > trace_change;
		trace_change = CalcTraceChange(tracers, cells, tess, fluxes, dt, hbc, interpolation, time, CellsEvolve,
			custom_evolution_manager, edge_velocities, lengths);
		MakeTracerExtensive(tracers, tess, cells, tracer_extensive);
		old_trace_extensive = tracer_extensive;
		UpdateTracerExtensive(tracer_extensive, trace_change, CellsEvolve, cells, tess, time);
	}

	vector<Conserved> intensive = CalcConservedIntensive(cells);

	vector<Conserved> extensive = CalcConservedExtensive(intensive, tess, pg);
	vector<Conserved> old_extensive = extensive;

	UpdateConservedExtensive(tess, fluxes, dt, extensive, hbc, lengths);

	ExternalForceContribution(tess, pg, cells, force, time, dt, extensive, hbc, fluxes, point_velocities, g, coldflows_flag, tracers, lengths);

	if (coldflows_flag)
	{
		Ek = GetMaxKineticEnergy(tess, cells, CellsEvolve);
		Ef = GetForceEnergy(tess, g);
	}

	vector<Primitive> old_cells = cells;
	vector<vector<double> > old_trace = tracers;

#ifndef RICH_MPI
	MoveMeshPoints(point_velocities, dt, tess);
#else
	if (procupdate != 0)
		procupdate->Update(proctess, tess);
	MoveMeshPoints(point_velocities, dt, tess, proctess);
#endif

#ifdef RICH_MPI
	vector<Conserved> ptoadd;
	vector<vector<double> > ttoadd;
	vector<size_t> ctemp(custom_evolution_indices);
	vector<size_t> ctoadd;
	SendRecvExtensive(extensive, tracer_extensive, custom_evolution_indices, tess.GetSentPoints(), tess.GetSentProcs(),
		ptoadd, ttoadd, ctoadd);

	KeepLocalPoints(extensive, tracer_extensive, custom_evolution_indices,
		tess.GetSelfPoint());

	if (!ptoadd.empty())
	{
		extensive.insert(extensive.end(), ptoadd.begin(),
			ptoadd.end());
	}

	if (!ttoadd.empty())
	{
		tracer_extensive.insert(tracer_extensive.end(), ttoadd.begin(), ttoadd.end());
	}

	if (!ctoadd.empty())
	{
		custom_evolution_indices.insert(custom_evolution_indices.end(), ctoadd.begin(),
			ctoadd.end());
	}

	SendRecvExtensive(old_extensive, old_trace_extensive, ctemp, tess.GetSentPoints(), tess.GetSentProcs(), ptoadd, ttoadd,
		ctoadd);

	KeepLocalPoints(old_extensive, old_trace_extensive, ctemp, tess.GetSelfPoint());

	if (!ptoadd.empty())
	{
		old_extensive.insert(old_extensive.end(), ptoadd.begin(),
			ptoadd.end());
	}

	if (!ttoadd.empty())
	{
		old_trace_extensive.insert(old_trace_extensive.end(), ttoadd.begin(), ttoadd.end());
	}

	vector<Vector2D> vel_toadd;
	SendRecvOldVector2D(point_velocities, tess.GetSentPoints(), tess.GetSentProcs(), vel_toadd);
	point_velocities = VectorValues(point_velocities, tess.GetSelfPoint());
	if (!vel_toadd.empty())
		point_velocities.insert(point_velocities.end(), vel_toadd.begin(), vel_toadd.end());

	CellsEvolve = convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	if (coldflows_flag)
	{
		vector<char> btoadd;
		SendRecvShockedStatus(shockedcells, tess.GetSentPoints(), tess.GetSentProcs(), btoadd);
		shockedcells = VectorValues(shockedcells, tess.GetSelfPoint());
		if (!btoadd.empty())
			shockedcells.insert(shockedcells.end(), btoadd.begin(), btoadd.end());
		vector<double> Ekadd;
		SendRecvVectorDouble(Ek, tess.GetSentPoints(), tess.GetSentProcs(),
			Ekadd);
		Ek = VectorValues(Ek, tess.GetSelfPoint());
		if (!Ekadd.empty())
			Ek.insert(Ek.end(), Ekadd.begin(), Ekadd.end());
		SendRecvVectorDouble(Ef, tess.GetSentPoints(), tess.GetSentProcs(),
			Ekadd);
		Ef = VectorValues(Ef, tess.GetSelfPoint());
		if (!Ekadd.empty())
			Ef.insert(Ef.end(), Ekadd.begin(), Ekadd.end());
	}
	if (densityfloor)
	{
		vector<Primitive> pmtoadd;
		ttoadd.clear();
		SendRecvPrimitive(old_cells, old_trace, tess.GetSentPoints(), tess.GetSentProcs(),
			pmtoadd, ttoadd, eos);
		KeepLocalPoints(old_cells, old_trace, tess.GetSelfPoint());
		if (!pmtoadd.empty())
			old_cells.insert(old_cells.end(), pmtoadd.begin(), pmtoadd.end());
		if (!ttoadd.empty())
			old_trace.insert(old_trace.begin(), ttoadd.begin(), ttoadd.end());
	}
#endif

	UpdateConservedIntensive(tess, extensive, intensive);

	if (coldflows_flag)
		FixPressure(intensive, tracer_extensive, eos, Ek, Ef, as, bs, CellsEvolve, tess,/*extensive,*/shockedcells, densityfloor);

	cells.resize(static_cast<size_t>(tess.GetPointNo()));
	vector<bool> min_density_on = UpdatePrimitives(intensive, eos, cells, CellsEvolve, old_cells, densityfloor,
		densitymin, pressuremin, tess, time + dt, tracers);
	if (traceflag)
		MakeTracerIntensive(tracers, tracer_extensive, tess, cells, pg, min_density_on, old_trace, CellsEvolve);

	// End first step

	// create ghost points if needed
#ifndef RICH_MPI
	PeriodicUpdateCells(cells, tracers, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetTotalPointNumber());

	PeriodicVelocityExchange(tess.GetAllCM(), tess.GetDuplicatedPoints(), tess.GetTotalPointNumber());
#else
	SendRecvHydro(cells, custom_evolution_indices, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), eos, tess.GetGhostIndeces(), tess.GetTotalPointNumber());

	SendRecvVelocity(point_velocities, tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());

	SendRecvVelocity(tess.GetAllCM(), tess.GetDuplicatedPoints(),
		tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
#endif

	CellsEvolve = convert_indices_to_custom_evolution(custom_evolution_manager,
		custom_evolution_indices);

	edge_velocities = tess.calc_edge_velocities(&hbc, point_velocities, time + dt);

	if (coldflows_flag&&EntropyCalc)
	{
		const int n = tess.GetPointNo();
		for (size_t i = 0; i<static_cast<size_t>(n); ++i)
			tracers[i][0] = eos.dp2s(cells[i].Density, cells[i].Pressure);
#ifdef RICH_MPI
		SendRecvTracers(tracers, tess.GetDuplicatedPoints(), tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
		if (tracers.size() != cells.size())
		{
			UniversalError eo("Tracers and cells have different length at second half time step");
			eo.AddEntry("CPU rank", get_mpi_rank());
			throw eo;
		}

#endif
	}

#ifdef RICH_MPI
	if (traceflag && (!coldflows_flag || !EntropyCalc))
	{
		SendRecvTracers(tracers, tess.GetDuplicatedPoints(),
			tess.GetDuplicatedProcs(), tess.GetGhostIndeces(), tess.GetTotalPointNumber());
		if (tracers.size() != cells.size())
		{
			UniversalError eo("Tracers and cells have different length at second half time step");
			eo.AddEntry("CPU rank", get_mpi_rank());
			throw eo;
		}
	}
#endif
	fluxes = calc_fluxes(tess, cells, dt, time + dt, interpolation, edge_velocities, hbc, rs, CellsEvolve,
		custom_evolution_manager, tracers);

	if (coldflows_flag)
	{
		const int n = tess.GetPointNo();
		shockedcells.resize(static_cast<size_t>(n));
		for (int i = 0; i<n; ++i)
			shockedcells[static_cast<size_t>(i)] = IsShockedCell(tess, i, cells, hbc, time) ? 1 : 0;
	}


	lengths = serial_generate(EdgeLengthCalculator(tess, pg));

	if (traceflag)
	{
		vector<vector<double> > trace_change;
		trace_change = CalcTraceChange(tracers, cells, tess, fluxes, dt, hbc, interpolation, time, CellsEvolve,
			custom_evolution_manager, edge_velocities, lengths);
		UpdateTracerExtensive(old_trace_extensive, trace_change, CellsEvolve, cells, tess, time);
	}

	UpdateConservedExtensive(tess, fluxes, dt, old_extensive, hbc, lengths);

	ExternalForceContribution(tess, pg, cells, force, time + dt, dt, old_extensive, hbc, fluxes, point_velocities, g, coldflows_flag, tracers, lengths);

	if (coldflows_flag)
	{
		Ek = GetMaxKineticEnergy(tess, cells, CellsEvolve);
		Ef = GetForceEnergy(tess, g);
	}

	extensive = 0.5*(old_extensive + extensive);
	tracer_extensive = 0.5*(tracer_extensive + old_trace_extensive);

	UpdateConservedIntensive(tess, extensive, intensive);

	if (coldflows_flag)
		FixPressure(intensive, old_trace_extensive, eos, Ek, Ef, as, bs, CellsEvolve, tess,
		/*extensive,*/shockedcells, densityfloor);

	min_density_on = UpdatePrimitives(intensive, eos, cells, CellsEvolve, old_cells, densityfloor,
		densitymin, pressuremin, tess, time + dt, tracers);

	if (traceflag)
	{
		MakeTracerIntensive(tracers, tracer_extensive, tess, cells, pg, min_density_on, old_trace, CellsEvolve);
	}
	return dt;
}
