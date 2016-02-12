#include "AreaFix.hpp"
#include <utility>
#include <set>

namespace
{
	Extensive ComputationalCell2Extensive(ComputationalCell const& cell, double volume, EquationOfState
		const& eos)
	{
		Extensive res;
		res.mass = cell.density*volume;
		res.energy = eos.dp2e(cell.density, cell.pressure, cell.tracers)*res.mass +
			0.5*res.mass*ScalarProd(cell.velocity, cell.velocity);
		res.momentum = res.mass*cell.velocity;
		for (size_t i = 0; i < cell.tracers.size(); ++i)
			res.tracers[i] = res.mass*cell.tracers[i];
		return res;
	}

	int GetEdgeIndex(Tessellation const& tessnew, int n0, int n1, int cell_index)
	{
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

	Vector2D FixPeriodicLeap(Vector2D const& point, Vector2D const& vel, double dt, OuterBoundary const& outer)
	{
		Vector2D res(0, 0);
		if ((point.x + vel.x*dt)>outer.GetGridBoundary(Right))
			res.x = outer.GetGridBoundary(Right) - outer.GetGridBoundary(Left);
		if ((point.x + vel.x*dt) < outer.GetGridBoundary(Left))
			res.x = outer.GetGridBoundary(Left) - outer.GetGridBoundary(Right);
		if ((point.y + vel.y*dt) > outer.GetGridBoundary(Up))
			res.y = outer.GetGridBoundary(Up) - outer.GetGridBoundary(Down);
		if ((point.y + vel.y*dt) < outer.GetGridBoundary(Down))
			res.y = outer.GetGridBoundary(Down) - outer.GetGridBoundary(Up);
		return res;
	}

	void GetAdjacentNeigh(Tessellation const& tessold, Edge const& edge, int cell_index,
		Vector2D const& added, int &n0other, int &n1other)
	{
		vector<int> const& edges = tessold.GetCellEdges(cell_index);
		vector<double> dist;
		Vector2D point;
		dist.reserve(edges.size());
		for (size_t i = 0; i < edges.size(); ++i)
		{
			Edge const& e = tessold.GetEdge(edges[i]);
			if (e.neighbors.first == edge.neighbors.first&&e.neighbors.second == edge.neighbors.second)
				dist.push_back(edge.GetLength() * 1000);
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

	int NewEdgeIndex(Tessellation const& tessold, Tessellation const& tessnew, int cell_index,
		Edge const& edge, int other_index)
	{
		int n0, n1;
		int res = -1;
		bool rigid = tessold.GetOriginalIndex(edge.neighbors.first) == tessold.GetOriginalIndex(edge.neighbors.second) ?
			true : false;
		// Find n0 and n1 the new flip neighbors
		GetAdjacentNeigh(tessold, edge, cell_index, Vector2D(0, 0), n0, n1);
		// Do we have two neighboring rigid edges?
		if (rigid && (n0 == cell_index || n1 == cell_index))
			return -1;

		// Are we flipping a real edge to a rigid one?
		if (n0 == cell_index)
		{
			n0 = n1;
			n1 = cell_index;
		}
		if (n1 == cell_index)
		{
			vector<int> const& newedges = tessnew.GetCellEdges(n0);
			for (size_t i = 0; i < newedges.size(); ++i)
			{
				Edge const& newedge = tessnew.GetEdge(newedges[i]);
				if (tessnew.GetOriginalIndex(newedge.neighbors.first) ==
					tessnew.GetOriginalIndex(newedge.neighbors.second))
					return newedges[i];
			}
			return -1;
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
		if (rigid)
			return res;
		int n0new, n1new;
		GetAdjacentNeigh(tessnew, tessnew.GetEdge(res), n0, Vector2D(0, 0), n0new, n1new);
		if ((n0new == cell_index&&n1new == other_index) || (n1new == cell_index&&other_index == n0new))
			return res;
		// We didn't find a new edge
		return -1;
	}

	bool FirstVertice(Edge edgenew, Tessellation const& tessnew, int cell, Vector2D const& edge_added)
	{
		edgenew.vertices.first += edge_added;
		edgenew.vertices.second += edge_added;
		const double R = 1e-8*tessnew.GetWidth(cell);
		vector<int> const& edges = tessnew.GetCellEdges(cell);
		for (size_t i = 0; i < edges.size(); ++i)
		{
			Edge edge = tessnew.GetEdge(edges[i]);
			if (edge.vertices.first.distance(edgenew.vertices.first) < R)
				return true;
			if (edge.vertices.first.distance(edgenew.vertices.second) < R)
				return false;
			if (edge.vertices.second.distance(edgenew.vertices.first) < R)
				return true;
			if (edge.vertices.second.distance(edgenew.vertices.second) < R)
				return false;
		}
		throw UniversalError("Couldn't find first vertice");
	}

	Edge FindOtherEdge(Tessellation const& tess, int cell_index, int other)
	{
		vector<int> const& edges = tess.GetCellEdges(cell_index);
		for (size_t i = 0; i < edges.size(); ++i)
		{
			Edge const& edge = tess.GetEdge(edges[i]);
			if (tess.GetOriginalIndex(edge.neighbors.first) == other)
				return edge;
			if (tess.GetOriginalIndex(edge.neighbors.second) == other)
				return edge;
		}
		throw UniversalError("Couldn't find other edge");
	}

	std::pair<Vector2D, Vector2D> TrianglesArea(Edge const& eold, Edge const &enew, Tessellation const& tessnew,
		Tessellation const& tessold, Vector2D const& new_edge_added, int cell_index, Vector2D const& cell_index_added)
	{ // first is the change in the old and the second the change in the new
	  // x is e1.first y is e1.second
		std::pair<Vector2D, Vector2D> res;
		boost::array<Vector2D, 4> points;
		double area_scale = 0;
		double sum = 0;
		bool rigid = tessold.GetOriginalIndex(eold.neighbors.first) ==
			tessold.GetOriginalIndex(eold.neighbors.second);
		Vector2D oldmeshpoint = (tessold.GetPointNo() < eold.neighbors.first) && rigid ?
			tessold.GetMeshPoint(eold.neighbors.second) : tessold.GetMeshPoint(eold.neighbors.first);
		TripleConstRef<Vector2D> temp(eold.vertices.first, eold.vertices.second, oldmeshpoint);
		if (orient2d(temp) > 0)
		{
			points[0] = eold.vertices.first;
			points[1] = eold.vertices.second;
		}
		else
		{
			points[1] = eold.vertices.first;
			points[0] = eold.vertices.second;
		}
		Vector2D  newmeshpoint = tessnew.GetMeshPoint(enew.neighbors.first);
		TripleConstRef<Vector2D> temp2(enew.vertices.first, enew.vertices.second, newmeshpoint);
		if (orient2d(temp2) > 0)
		{
			points[3] = enew.vertices.first + new_edge_added;
			points[2] = enew.vertices.second + new_edge_added;
		}
		else
		{
			points[2] = enew.vertices.first + new_edge_added;
			points[3] = enew.vertices.second + new_edge_added;
		}
		int npoints = tessold.GetPointNo();
		Vector2D toadd = new_edge_added - cell_index_added;
		bool first = FirstVertice(enew, tessnew, cell_index, toadd);
		if (eold.neighbors.first > npoints)
		{
			first = !first;
			if (rigid)
			{
				Vector2D tempp = points[0];
				points[0] = points[1];
				points[1] = tempp;
			}
		}
		Vector2D third_point = (first ? enew.vertices.first : enew.vertices.second) + new_edge_added;
		{
			TripleConstRef<Vector2D> temp3 = TripleConstRef<Vector2D>(points[0], third_point, points[1]);
			res.first.x = 0.5*orient2d(temp3);
			area_scale = std::abs(res.first.x);
			sum = res.first.x;
		}
		{
			third_point = (first ? enew.vertices.second : enew.vertices.first) + new_edge_added;
			TripleConstRef<Vector2D> temp4 = TripleConstRef<Vector2D>(points[1], third_point, points[0]);
			res.first.y = 0.5*orient2d(temp4);
			area_scale = std::max(area_scale, std::abs(res.first.y));
			sum += res.first.y;
		}
		bool newrigid = tessnew.GetOriginalIndex(enew.neighbors.first) == tessnew.GetOriginalIndex(enew.neighbors.second);
		int other = tessnew.GetOriginalIndex(enew.neighbors.first);
		Edge otheredge = FindOtherEdge(tessold, cell_index, other);
		Vector2D added;
		if (tessold.GetOriginalIndex(otheredge.neighbors.first) == other)
			added = tessold.GetMeshPoint(other) - tessold.GetMeshPoint(otheredge.neighbors.first);
		else
			added = tessold.GetMeshPoint(other) - tessold.GetMeshPoint(otheredge.neighbors.second);
		first = FirstVertice(eold, tessold, other, added);

		if (newrigid&&enew.neighbors.first > npoints)
			first = !first;
		if (enew.neighbors.first < npoints || (tessnew.GetOriginalIndex(enew.neighbors.first) !=
			tessnew.GetOriginalIndex(enew.neighbors.second)))
		{
			third_point = first ? eold.vertices.first : eold.vertices.second;
			TripleConstRef<Vector2D> temp5(points[2], third_point, points[3]);
			res.second.x = 0.5*orient2d(temp5);
			area_scale = std::max(area_scale, std::abs(res.second.x));
			sum += res.second.x;
		}

		if (enew.neighbors.second < npoints || (tessnew.GetOriginalIndex(enew.neighbors.first) !=
			tessnew.GetOriginalIndex(enew.neighbors.second)))
		{
			third_point = !first ? eold.vertices.first : eold.vertices.second;
			TripleConstRef<Vector2D> temp6(points[3], third_point, points[2]);
			res.second.y = 0.5*orient2d(temp6);
			area_scale = std::max(area_scale, std::abs(res.second.y));
			sum += res.second.y;
		}
		if (std::abs(sum) > 1e-5*area_scale)
		{
			UniversalError eo("Bad total area in TrianglesArea");
			eo.AddEntry("area_scale", area_scale);
			eo.AddEntry("sum", sum);
			eo.AddEntry("old neigh 0", eold.neighbors.first);
			eo.AddEntry("old neigh 1", eold.neighbors.second);
			throw eo;
		}
		return res;
	}

	Vector2D FixNewEdgeLeap(Tessellation const& tessold, int real_cell, int cell_index)
	{
		vector<int> const& edges = tessold.GetCellEdges(cell_index);
		for (size_t i = 0; i < edges.size(); ++i)
		{
			Edge const& edge = tessold.GetEdge(edges[i]);
			if (tessold.GetOriginalIndex(edge.neighbors.first) == real_cell)
				return tessold.GetMeshPoint(edge.neighbors.first) - tessold.GetMeshPoint(real_cell);
			if (tessold.GetOriginalIndex(edge.neighbors.second) == real_cell)
				return tessold.GetMeshPoint(edge.neighbors.second) - tessold.GetMeshPoint(real_cell);
		}
		throw UniversalError("Couldn't find New edge leap");
	}

	void GetCellIndex(Edge const& edge, Tessellation const& tess, int &cell_index, int &other_index)
	{
		int npoints = tess.GetPointNo();
		if (edge.neighbors.first < npoints)
		{
			cell_index = edge.neighbors.first;
			other_index = tess.GetOriginalIndex(edge.neighbors.second);
		}
		else
		{
			cell_index = edge.neighbors.second;
			other_index = tess.GetOriginalIndex(edge.neighbors.first);
		}
	}

	double GetdAFlux(int mid_index, Tessellation const& tessmid, double dt,
		int other_index, vector<Vector2D> const& fv)
	{
		double res = 0;
		if (mid_index >= 0)
		{
			Edge const& edge2 = tessmid.GetEdge(mid_index);
			const Vector2D norm = tessmid.GetMeshPoint(edge2.neighbors.second) - tessmid.GetMeshPoint(edge2.neighbors.first);
			res = ScalarProd(norm, fv[static_cast<size_t>(mid_index)])*dt*edge2.GetLength() / abs(norm);
			if (tessmid.GetOriginalIndex(edge2.neighbors.first) == other_index)
				res *= -1;
		}
		return res;
	}

	void AreaFixEdgeDisappear(Tessellation const& tessold, Tessellation const& tessnew, int cell_index, int other_index,
		Edge const& edge, OuterBoundary const& outer, int npoints, vector<Vector2D> const& pointvelocity, double dt,
		double dA_flux, int mid_index, Tessellation const& tessmid, vector<Vector2D> const& fv,
		vector<ComputationalCell> const& cells, EquationOfState const& eos, vector<Extensive> &res,
		std::set<std::pair<int, int > > &flipped_set, Vector2D const& cell_added)
	{
		// Is there a corresponding new edge? An edge flip
		int new_edge = NewEdgeIndex(tessold, tessnew, cell_index, edge, other_index);
		if (new_edge < 0)
			return;
		bool both_real = edge.neighbors.first < npoints && edge.neighbors.second < npoints;
		std::pair<int, int> neighbors(tessold.GetOriginalIndex(edge.neighbors.first), tessold.GetOriginalIndex(edge.neighbors.second));
		if (flipped_set.count(neighbors) > 0)
			return;
		Edge edge_new = tessnew.GetEdge(new_edge);
		// Did it jump?
		Vector2D addedother;
		if (outer.GetBoundaryType() == Periodic)
		{
			if (edge_new.neighbors.first < npoints)
			{
				addedother = -1 * FixPeriodicLeap(tessnew.GetMeshPoint(edge_new.neighbors.first),
					pointvelocity[static_cast<size_t>(edge_new.neighbors.first)], -dt, outer);
				addedother += FixNewEdgeLeap(tessold, edge_new.neighbors.first, cell_index);
			}
			else
			{
				addedother = -1 * FixPeriodicLeap(tessnew.GetMeshPoint(edge_new.neighbors.second),
					pointvelocity[static_cast<size_t>(edge_new.neighbors.second)], -dt, outer);
				addedother += FixNewEdgeLeap(tessold, edge_new.neighbors.second, cell_index);
			}
		}
		std::pair<Vector2D, Vector2D> areas = TrianglesArea(edge, edge_new, tessnew, tessold, addedother, cell_index,
			cell_added);
		vector<double> dA(4);
		if (mid_index >= 0 && tessmid.GetOriginalIndex(tessmid.GetEdge(mid_index).neighbors.first) == other_index)
			dA_flux *= -1;
		if (mid_index >= 0 && tessmid.GetOriginalIndex(tessmid.GetEdge(mid_index).neighbors.first) == tessold.GetOriginalIndex(edge.neighbors.second))
			dA_flux *= -1;
		dA[0] = areas.first.x - dA_flux;
		dA[1] = areas.first.y + dA_flux;
		dA[2] = areas.second.x;
		dA[3] = areas.second.y;

		// Did the new edge exist in the mid time step? This should mean mid_index=-1
		int mid_index2 = GetEdgeIndex(tessmid, tessnew.GetOriginalIndex(edge_new.neighbors.first),
			tessnew.GetOriginalIndex(edge_new.neighbors.second), tessnew.GetOriginalIndex(edge_new.neighbors.first));

		if (tessnew.GetOriginalIndex(edge_new.neighbors.first) == tessnew.GetOriginalIndex(edge_new.neighbors.second))
			mid_index2 = -1;

		double dA_flux2 = GetdAFlux(mid_index2, tessmid, dt, tessnew.GetOriginalIndex(edge_new.neighbors.first), fv);
		dA[2] += dA_flux2;
		dA[3] -= dA_flux2;

		// The conserved to remove
		Extensive TotalRemoved(cells[0].tracers);
		double total_added_area = 0;
		if (dA[0] < 0)
		{
			Extensive toremove = ComputationalCell2Extensive(cells[
				static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.first))], -dA[0], eos);
			res[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.first))] -= toremove;
			TotalRemoved += toremove;
		}
		else
			total_added_area += dA[0];
		if (dA[1] < 0)
		{
			Extensive toremove = ComputationalCell2Extensive(cells[
				static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.second))], -dA[1], eos);
			res[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.second))] -= toremove;
			TotalRemoved += toremove;
		}
		else
			total_added_area += dA[1];

		if (dA[2] < 0)
		{
			Extensive toremove = ComputationalCell2Extensive(cells[
				static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))], -dA[2], eos);
			res[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))] -= toremove;
			TotalRemoved += toremove;
		}
		else
			total_added_area += dA[2];
		if (dA[3] < 0)
		{
			Extensive toremove = ComputationalCell2Extensive(cells[
				static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))], -dA[3], eos);
			res[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))] -= toremove;
			TotalRemoved += toremove;
		}
		else
			total_added_area += dA[3];

		// Add the conserved
		if (dA[0] > 0)
			res[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.first))] += (dA[0] / total_added_area)*TotalRemoved;
		if (dA[1] > 0)
			res[static_cast<size_t>(tessold.GetOriginalIndex(edge.neighbors.second))] += (dA[1] / total_added_area)*TotalRemoved;
		if (dA[2] > 0)
			if ((edge_new.neighbors.first < tessnew.GetPointNo()) ||
				(tessnew.GetOriginalIndex(edge_new.neighbors.first) != tessnew.GetOriginalIndex(edge_new.neighbors.second)))
				res[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.first))] += (dA[2] / total_added_area)*TotalRemoved;
		if (dA[3] > 0)
			if ((edge_new.neighbors.second < tessnew.GetPointNo()) ||
				(tessnew.GetOriginalIndex(edge_new.neighbors.first) != tessnew.GetOriginalIndex(edge_new.neighbors.second)))
				res[static_cast<size_t>(tessnew.GetOriginalIndex(edge_new.neighbors.second))] += (dA[3] / total_added_area)*TotalRemoved;
		if (!both_real)
		{
			flipped_set.insert(neighbors);
			std::pair<int, int> neigh_temp(neighbors.second, neighbors.first);
			flipped_set.insert(neigh_temp);
		}
	}

	double EdgesArea(Edge const& e1, Edge const &e2, Vector2D const& pointold, Vector2D const& pointnew)
	{
		// e1 is old tess, e2 is new tess
		boost::array<Vector2D, 4> edge_points;
		double res = 0;
		TripleConstRef<Vector2D> temp1(pointold, e1.vertices.first, e1.vertices.second);
		TripleConstRef<Vector2D> temp2(pointnew, e2.vertices.first, e2.vertices.second);
		if (orient2d(temp1) > 0)
		{
			edge_points[0] = e1.vertices.first;
			edge_points[1] = e1.vertices.second;
		}
		else
		{
			edge_points[1] = e1.vertices.first;
			edge_points[0] = e1.vertices.second;
		}
		if (orient2d(temp2) > 0)
		{
			edge_points[2] = e2.vertices.first;
			edge_points[3] = e2.vertices.second;
		}
		else
		{
			edge_points[3] = e2.vertices.first;
			edge_points[2] = e2.vertices.second;
		}
		TripleConstRef<Vector2D> tri1(edge_points[0], edge_points[3], edge_points[1]);
		TripleConstRef<Vector2D> tri2(edge_points[0], edge_points[2], edge_points[3]);
		res += orient2d(tri1);
		res += orient2d(tri2);
		return res*0.5;
	}

	ComputationalCell GetDonorCell(ComputationalCell const& c0, ComputationalCell const& c1,
		bool seconddonor, EquationOfState const& /*eos*/)
	{
		ComputationalCell p_mid;
		double ratio = c1.density / c0.density;
		if (ratio<2 && ratio>0.5)
			p_mid.density = 0.5*(c1.density + c0.density);
		else
		{
			if (seconddonor)
				p_mid.density = c1.density;
			else
				p_mid.density = c0.density;
		}
		ratio = c1.pressure / c0.pressure;
		if (ratio<2 && ratio>0.5)
			p_mid.pressure = 0.5*(c1.pressure + c0.pressure);
		else
		{
			if (seconddonor)
				p_mid.pressure = c1.pressure;
			else
				p_mid.pressure = c0.pressure;
		}
		ratio = c1.velocity.x / c0.velocity.x;
		if (ratio<2 && ratio>0.5)
			p_mid.velocity.x = 0.5*(c1.velocity.x + c0.velocity.x);
		else
		{
			if (seconddonor)
				p_mid.velocity.x = c1.velocity.x;
			else
				p_mid.velocity.x = c0.velocity.x;
		}
		ratio = c1.velocity.y / c0.velocity.y;
		if (ratio<2 && ratio>0.5)
			p_mid.velocity.y = 0.5*(c1.velocity.y + c0.velocity.y);
		else
		{
			if (seconddonor)
				p_mid.velocity.y = c1.velocity.y;
			else
				p_mid.velocity.y = c0.velocity.y;
		}
		for (size_t i = 0; i < c0.tracers.size(); ++i)
		{
			ratio = c0.tracers[i] / c1.tracers[i];
			if (ratio<2 && ratio>0.5)
				p_mid.tracers[i] = 0.5*c0.tracers[i] * (1 + 1 / ratio);
			else
				if (seconddonor)
					p_mid.tracers[i] = c1.tracers[i];
				else
					p_mid.tracers[i] = c0.tracers[i];
		}
		return p_mid;
	}
}

vector<Extensive> FluxFix2(Tessellation const& tessold, Tessellation const& tessmid,
	Tessellation const& tessnew, vector<Vector2D> const& pointvelocity, double dt,
	vector<ComputationalCell> const& cells, vector<Extensive> const& /*fluxes*/,
	vector<Vector2D> const& fv, OuterBoundary const& outer, EquationOfState const& eos)
{
	size_t npoints = static_cast<size_t>(tessold.GetPointNo());
	vector<Extensive> res(static_cast<size_t>(npoints), Extensive(cells[0].tracers));
	// Fix the fluxes
	Vector2D temp0, temp1;
	std::set<std::pair<int, int > > flipped_set, only_mid;
	size_t nedgesold = static_cast<size_t>(tessold.GetTotalSidesNumber());
	for (size_t i = 0; i < nedgesold; ++i)
	{
		Edge const& edge = tessold.GetEdge(static_cast<int>(i));

		int cell_index, other_index;
		GetCellIndex(edge, tessold, cell_index, other_index);
		int mid_index = GetEdgeIndex(tessmid, cell_index, other_index, cell_index);
		int new_index = GetEdgeIndex(tessnew, cell_index, other_index, cell_index);
		double dA_flux = GetdAFlux(mid_index, tessmid, dt, other_index, fv);

		Vector2D real_p;
		real_p = tessold.GetMeshPoint(cell_index);
		Vector2D added(0, 0);
		// check if periodic jumped
		if (outer.GetBoundaryType() == Periodic)
			added = FixPeriodicLeap(real_p, pointvelocity[static_cast<size_t>(cell_index)], dt, outer);

		// Did the edge disappear? An edge flip
		if (new_index < 0)
		{
			AreaFixEdgeDisappear(tessold, tessnew, cell_index, other_index, edge, outer, static_cast<int>(npoints), pointvelocity, dt,
				dA_flux, mid_index, tessmid, fv, cells, eos, res, flipped_set, added);
			continue;
		}
		Vector2D real_p_new = tessnew.GetMeshPoint(cell_index);
		Edge edge_other = tessnew.GetEdge(new_index);
		if (outer.GetBoundaryType() == Periodic)
		{
			edge_other.vertices.first += added;
			edge_other.vertices.second += added;
			real_p_new += added;
		}
		double area = EdgesArea(edge, edge_other, real_p, real_p_new);
		ComputationalCell p_mid = GetDonorCell(cells[static_cast<size_t>(cell_index)], cells[static_cast<size_t>(other_index)],
			(area - dA_flux)>0, eos);
		Extensive toadd = ComputationalCell2Extensive(p_mid, (-dA_flux + area), eos);

		if (other_index == edge.neighbors.first || other_index == edge.neighbors.second)
			res[static_cast<size_t>(other_index)] -= toadd;
		if (cell_index == edge.neighbors.first || cell_index == edge.neighbors.second)
			res[static_cast<size_t>(cell_index)] += toadd;
		if (mid_index < 0)
		{
			// Was there a flip in mid step?
			int mid_edge = NewEdgeIndex(tessold, tessmid, cell_index, edge, other_index);
			if (mid_edge < 0)
				continue;
			Edge const& edge2 = tessmid.GetEdge(mid_edge);
			if (tessmid.GetOriginalIndex(edge2.neighbors.first) == tessmid.GetOriginalIndex(edge2.neighbors.second))
				continue;
			double dA_flux2 = GetdAFlux(mid_edge, tessmid, dt, tessmid.GetOriginalIndex(edge2.neighbors.second), fv);
			Extensive toadd2(cells[0].tracers);
			if (dA_flux2 > 0)
				toadd2 = ComputationalCell2Extensive(cells[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.second))],
					dA_flux2, eos);
			else
				toadd2 = ComputationalCell2Extensive(cells[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.first))],
					dA_flux2, eos);
			std::pair<int, int> mid_neigh(tessmid.GetOriginalIndex(edge2.neighbors.first), tessmid.GetOriginalIndex(edge2.neighbors.second));
			std::pair<int, int> mid_neigh2(mid_neigh.second, mid_neigh.first);
			if (only_mid.count(mid_neigh) == 0)
			{
				only_mid.insert(mid_neigh);
				only_mid.insert(mid_neigh2);
				res[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.second))] += toadd2;
				res[static_cast<size_t>(tessmid.GetOriginalIndex(edge2.neighbors.first))] -= toadd2;
			}
		}
	}
	return res;
}
