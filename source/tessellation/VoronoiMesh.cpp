#include "VoronoiMesh.hpp"
#include <cmath>
#include "../misc/simple_io.hpp"
#include "hdf5_logger.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#include "../mpi/mpi_commands.hpp"
#endif

using std::abs;
using std::endl;
using std::cout;

namespace
{
#ifdef RICH_MPI
	template<class T> void tidy(vector<T>& v)
	{
		if (!v.empty()) {
			sort(v.begin(), v.end());
			v = unique(v);
		}
	}
#endif

	vector<Vector2D> UpdatePoints(vector<Vector2D> const& points, OuterBoundary const* obc)
	{
		if (obc->GetBoundaryType() == Rectengular)
			return points;
		vector<Vector2D> res;
		res.reserve(points.size());
		int npoints = static_cast<int>(points.size());
		const double dx = obc->GetGridBoundary(Right) - obc->GetGridBoundary(Left);
		const double dy = obc->GetGridBoundary(Up) - obc->GetGridBoundary(Down);
		for (int i = 0; i < npoints; ++i)
		{
			Vector2D temp(points[static_cast<size_t>(i)]);
			if (obc->GetBoundaryType() == Periodic)
			{
				if (temp.x > obc->GetGridBoundary(Right))
					temp.x -= dx;
				if (temp.x < obc->GetGridBoundary(Left))
					temp.x += dx;
				if (temp.y > obc->GetGridBoundary(Up))
					temp.y -= dy;
				if (temp.y < obc->GetGridBoundary(Down))
					temp.y += dy;
			}
			if (obc->GetBoundaryType() == HalfPeriodic)
			{
				if (temp.x > obc->GetGridBoundary(Right))
					temp.x -= dx;
				if (temp.x < obc->GetGridBoundary(Left))
					temp.x += dx;
			}
			res.push_back(temp);
		}
		return res;
	}
}

VoronoiMesh::VoronoiMesh
(vector<Vector2D> const& points,
	OuterBoundary const& bc, bool HOrder) :
	logger(0),
	eps(1e-8),
	obc(0),
	cell_edges(vector<Edge>()),
	edges(vector<Edge>()),
	CM(vector<Vector2D>()),
	mesh_vertices(vector<vector<int> >()),
	Tri(),
	GhostProcs(vector<int>()),
	GhostPoints(vector<vector<int> >()),
	SentProcs(vector<int>()),
	SentPoints(vector<vector<int> >()),
	selfindex(vector<size_t>()),
	NGhostReceived(vector<vector<int> >()),
	OrgCorner(),
	Nextra(0)
{
	Initialise(points, &bc, HOrder);
}

#ifdef RICH_MPI
VoronoiMesh::VoronoiMesh
(Tessellation const& proctess,
	vector<Vector2D> const& points,
	OuterBoundary const& bc,
	bool HOrder) :
	logger(0),
	eps(1e-8),
	obc(0),
	cell_edges(vector<Edge>()),
	edges(vector<Edge>()),
	CM(vector<Vector2D>()),
	mesh_vertices(vector<vector<int> >()),
	Tri(),
	GhostProcs(vector<int>()),
	GhostPoints(vector<vector<int> >()),
	SentProcs(vector<int>()),
	SentPoints(vector<vector<int> >()),
	selfindex(vector<size_t>()),
	NGhostReceived(vector<vector<int> >()),
	OrgCorner(),
	Nextra(0)
{
	Initialise(points, proctess, &bc, HOrder);
}
#endif

vector<int> VoronoiMesh::AddPointsAlongEdge(size_t point, vector<vector<int> > const&copied,
	int side)
{
	int ncopy = static_cast<int>(copied[static_cast<size_t>(side)].size());
	Vector2D vec = Tri.get_point(point);
	vector<double> dist(static_cast<size_t>(ncopy));
	for (size_t i = 0; i < copied[static_cast<size_t>(side)].size(); ++i)
		dist[i] = vec.distance(Tri.get_point(static_cast<size_t>(copied[static_cast<size_t>(side)][i])));
	const int copylength = min(7, static_cast<int>(copied[static_cast<size_t>(side)].size()) - 1);
	vector<int> index, toadd(static_cast<size_t>(copylength));
	sort_index(dist, index);
	for (int i = 0; i < copylength; ++i)
		toadd[static_cast<size_t>(i)] = copied[static_cast<size_t>(side)][static_cast<size_t>(index[static_cast<size_t>(i) + 1])];
	return toadd;
}

Vector2D VoronoiMesh::CalcFaceVelocity(Vector2D wl, Vector2D wr, Vector2D rL, Vector2D rR,
	Vector2D f)const
{
	const Vector2D wprime = ScalarProd(wl - wr, f - (rR + rL) / 2)*(rR - rL) / pow(abs(rR - rL), 2);
	return 0.5*(wl + wr) + wprime;
}

bool VoronoiMesh::NearBoundary(int index) const
{
	const int n = int(mesh_vertices[static_cast<size_t>(index)].size());
	const int N = Tri.get_length();
	for (int i = 0; i < n; ++i)
	{
		const int n0 = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].neighbors.first;
		const int n1 = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].neighbors.second;
		if (n0 < 0 || n1 < 0 || n0 >= N || n1 >= N)
			return true;
	}
	return false;
}

VoronoiMesh::~VoronoiMesh(void) {}

int VoronoiMesh::GetOriginalIndex(int point) const
{
	int npoints = GetPointNo();
	if (point < npoints)
		return point;
	else
		return Tri.GetOrgIndex(point);
}

vector<size_t> VoronoiMesh::GetSelfPoint(void)const
{
	return selfindex;
}

VoronoiMesh::VoronoiMesh(void) :
	logger(0),
	eps(1e-8),
	obc(0),
	cell_edges(vector<Edge>()),
	edges(vector<Edge>()),
	CM(vector<Vector2D>()),
	mesh_vertices(vector<vector<int> >()),
	Tri(),
	GhostProcs(vector<int>()),
	GhostPoints(vector<vector<int> >()),
	SentProcs(vector<int>()),
	SentPoints(vector<vector<int> >()),
	selfindex(vector<size_t>()),
	NGhostReceived(vector<vector<int> >()),
	OrgCorner(),
	Nextra(0) {}

VoronoiMesh::VoronoiMesh(VoronoiMesh const& other) :
	logger(other.logger),
	eps(other.eps),
	obc(other.obc),
	cell_edges(other.cell_edges),
	edges(other.edges),
	CM(other.CM),
	mesh_vertices(other.mesh_vertices),
	Tri(other.Tri),
	GhostProcs(other.GhostProcs),
	GhostPoints(other.GhostPoints),
	SentProcs(other.SentProcs),
	SentPoints(other.SentPoints),
	selfindex(other.selfindex),
	NGhostReceived(other.NGhostReceived),
	OrgCorner(),
	Nextra(other.Nextra)
{}

void VoronoiMesh::build_v()
{
	Vector2D center, center_temp;
	int j;
	facet to_check;
	Edge edge_temp;
	Vector2D p_temp;
	mesh_vertices.clear();
	mesh_vertices.resize(static_cast<size_t>(Tri.get_length()));
	edges.reserve(static_cast<size_t>(Tri.get_length()*3.5));
	int N = Tri.GetOriginalLength();
	for (int i = 0; i < N; ++i)
		mesh_vertices[static_cast<size_t>(i)].reserve(7);
	int Nfacets = Tri.get_num_facet();
	vector<Vector2D> centers(static_cast<size_t>(Nfacets));
	double R = 0;
	for (int i = 0; i < Nfacets; ++i)
		centers[static_cast<size_t>(i)] = Tri.GetCircleCenter(i,R);
	for (int i = 0; i < Nfacets; ++i)
	{
		center = centers[static_cast<size_t>(i)];
		to_check = Tri.get_facet(i);
		for (j = 0; j < 3; ++j)
		{
			if (to_check.neighbors[static_cast<size_t>(j)] == Tri.get_last_loc())
				continue;
			if (to_check.neighbors[static_cast<size_t>(j)] < i)
				continue;
			center_temp = centers[static_cast<size_t>(to_check.neighbors[static_cast<size_t>(j)])];
			{
				edge_temp.vertices.first = center;
				edge_temp.vertices.second = center_temp;
				edge_temp.neighbors.first = to_check.vertices[static_cast<size_t>(j)];
				edge_temp.neighbors.second = to_check.vertices[static_cast<size_t>(j + 1) % 3];
				if (legal_edge(&edge_temp))
				{
					// I added a change here, if edge has zero length I don't add it.
					if (edge_temp.GetLength() > eps*sqrt(Tri.GetFacetRadius(i)*
						Tri.GetFacetRadius(to_check.neighbors[static_cast<size_t>(j)])))
					{
						{
							if (edge_temp.neighbors.first < Tri.GetOriginalLength())
								mesh_vertices[static_cast<size_t>(edge_temp.neighbors.first)].push_back(static_cast<int>(edges.size()));
							if (edge_temp.neighbors.second < Tri.GetOriginalLength())
								mesh_vertices[static_cast<size_t>(edge_temp.neighbors.second)].push_back(static_cast<int>(edges.size()));
							edges.push_back(edge_temp);
						}
					}
				}
			}
		}
	}
}

namespace
{
	vector<Vector2D> calc_procpoints(const OuterBoundary& bc)
	{
		vector<Vector2D> res(4);
		res[0] = Vector2D(bc.GetGridBoundary(Left), bc.GetGridBoundary(Down));
		res[1] = Vector2D(bc.GetGridBoundary(Right), bc.GetGridBoundary(Down));
		res[2] = Vector2D(bc.GetGridBoundary(Right), bc.GetGridBoundary(Up));
		res[3] = Vector2D(bc.GetGridBoundary(Left), bc.GetGridBoundary(Up));
		return res;
	}

	/*
	Vector2D GetReflection(OuterBoundary const& bc, size_t index, Vector2D const& point)
	{
	switch (index)
	{
	case 0:
	return point + Vector2D(2 * (bc.GetGridBoundary(Right) - point.x), 0);
	case 1:
	return point + Vector2D(0, 2 * (bc.GetGridBoundary(Up) - point.y));
	case 2:
	return point + Vector2D(2 * (bc.GetGridBoundary(Left) - point.x), 0);
	case 3:
	return point + Vector2D(0, 2 * (bc.GetGridBoundary(Down) - point.y));
	}
	throw UniversalError("Wrong index in VoronoiMesh::GetReflection");
	}*/
}

void VoronoiMesh::Initialise(vector<Vector2D>const& pv, OuterBoundary const* _bc, bool reorder)
{
	obc = _bc;
	vector<Vector2D> points;
	if (reorder)
		points = VectorValues(pv, HilbertOrder(pv, static_cast<int>(pv.size())));
	else
		points = pv;
	Tri.build_delaunay(UpdatePoints(points, obc), calc_procpoints(*obc));

	Nextra = static_cast<int>(Tri.ChangeCor().size());
	Tri.BuildBoundary(_bc, _bc->GetBoxEdges());
	
	eps = 1e-8;
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	build_v();

	if (logger)
		logger->output(*this);

	CM.resize(Tri.getCor().size());
	for (size_t i = 0; i < pv.size(); ++i)
		CM[i] = CalcCellCM(i);

	size_t counter = pv.size() + 3;
	if (_bc->GetBoundaryType() == Periodic)
	{
		for (int i = static_cast<int>(counter); i < Tri.GetCorSize(); ++i)
		{
			int NorgIndex = Tri.GetOrgIndex(i);
			if (NorgIndex < Nextra)
			{
				CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + 
				  (Tri.get_point(static_cast<size_t>(i)) -
				   Tri.get_point(static_cast<size_t>(NorgIndex)));
			}
		}
	}
	else
	{
		if (_bc->GetBoundaryType() == Rectengular)
		{
			for (int i = static_cast<int>(counter); i < Tri.GetCorSize(); ++i)
			{
				int NorgIndex = Tri.GetOrgIndex(i);
				if (NorgIndex < Nextra)
				{
				  Vector2D norm = normalize(Tri.get_point(static_cast<size_t>(i)) - 
								  Tri.get_point(static_cast<size_t>(NorgIndex)));
				  CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + Tri.get_point(static_cast<size_t>(i))
						- Tri.get_point(static_cast<size_t>(NorgIndex)) - 2 * norm * ScalarProd(norm,
							CM[static_cast<size_t>(NorgIndex)] - Tri.get_point(static_cast<size_t>(NorgIndex)));
				}
			}
		}
		else  // Half periodic case
		{
			double dx = _bc->GetGridBoundary(Right) - _bc->GetGridBoundary(Left);
			for (int i = static_cast<int>(counter); i < Tri.GetCorSize(); ++i)
			{
				int NorgIndex = Tri.GetOrgIndex(i);
				if (NorgIndex < Nextra)
				{
				  double dx_temp = fastabs(Tri.get_point(static_cast<size_t>(i)) - Tri.get_point(static_cast<size_t>(NorgIndex)));
					if (dx_temp<1.0001*dx && dx_temp*1.0001>dx)
					  CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + (Tri.get_point(static_cast<size_t>(i)) -
							Tri.get_point(static_cast<size_t>(NorgIndex)));
					else
					{
					  Vector2D norm = normalize(Tri.get_point(static_cast<size_t>(i)) - Tri.get_point(static_cast<size_t>(NorgIndex)));
					  CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + Tri.get_point(static_cast<size_t>(i))
							- Tri.get_point(static_cast<size_t>(NorgIndex)) - 2 * norm * ScalarProd(norm,
								CM[static_cast<size_t>(NorgIndex)] - Tri.get_point(static_cast<size_t>(NorgIndex)));
					}
				}
			}
		}
	}
}

bool VoronoiMesh::legal_edge(Edge *e) //checks if both ends of the edge are outside the grid and that the edge doesn't cross the grid
{
	bool n0 = e->neighbors.first < Tri.get_length();
	bool n1 = e->neighbors.second < Tri.get_length();
	if ( n0 || n1)
	{
		if (obc->GetBoundaryType() == Rectengular)
		{
			if (!n1 || !n0)
			{
#ifdef RICH_MPI
				if (!n0)
				{
					Vector2D p = Tri.get_point(n0);
					if (p.x > obc->GetGridBoundary(Right) || p.x<obc->GetGridBoundary(Left) ||
						p.y>obc->GetGridBoundary(Up) || p.y < obc->GetGridBoundary(Down))
					{
						if (GetOriginalIndex(n0) != GetOriginalIndex(n1))
							return false;
					}
				}
				else
				{
					Vector2D p = Tri.get_point(n1);
					if (p.x > obc->GetGridBoundary(Right) || p.x<obc->GetGridBoundary(Left) ||
						p.y>obc->GetGridBoundary(Up) || p.y < obc->GetGridBoundary(Down))
					{
						if (GetOriginalIndex(n0) != GetOriginalIndex(n1))
							return false;
					}
				}
#else
				if (GetOriginalIndex(e->neighbors.first) != GetOriginalIndex(e->neighbors.second))
					return false;
#endif
			}
		}
		return true;
	}
	else
		return false;
}

double VoronoiMesh::GetWidth(int index)const
{
	return sqrt(GetVolume(index) / M_PI);
}

vector<int> const& VoronoiMesh::GetCellEdges(int index) const
{
	return mesh_vertices.at(static_cast<size_t>(index));
}

double VoronoiMesh::GetVolume(int index) const
{
	const Vector2D center = Tri.get_point(static_cast<size_t>(index));
	double area = 0;
	const size_t nloop = mesh_vertices[static_cast<size_t>(index)].size();
	for (size_t i = 0; i < nloop; ++i)
	{
		const Vector2D p1 = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].vertices.first -
			center;
		const Vector2D p2 = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].vertices.second -
			center;
		area += 0.5*abs(ScalarProd(p1, zcross(p2)));
	}
	return area;
}

Vector2D VoronoiMesh::CalcCellCM(size_t index) const
{
	const Vector2D center = edges[static_cast<size_t>(mesh_vertices[index].front())].vertices.first;
	Vector2D pc(0, 0);
	double area = 0;
	for (size_t i = 1; i < mesh_vertices[index].size(); i++)
	{
		const Edge& edge = edges[static_cast<size_t>(mesh_vertices[index][i])];
		const Vector2D p1 = edge.vertices.first - center;
		const Vector2D p2 = edge.vertices.second - center;
		const double area_temp = 0.5*abs(ScalarProd(p1, zcross(p2)));
		area += area_temp;
		pc += (area_temp / 3.)*(center + edge.vertices.first + edge.vertices.second);
	}
	return pc / area;
}

vector<Vector2D>& VoronoiMesh::GetMeshPoints(void)
{
	return Tri.GetMeshPoints();
}

vector<int> VoronoiMesh::Update(const vector<Vector2D>& pv, bool reorder)
{
	// Clean_up last step
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	vector<Vector2D> procpoints;
	procpoints.push_back(Vector2D(obc->GetGridBoundary(Left), obc->GetGridBoundary(Down)));
	procpoints.push_back(Vector2D(obc->GetGridBoundary(Right), obc->GetGridBoundary(Down)));
	procpoints.push_back(Vector2D(obc->GetGridBoundary(Right), obc->GetGridBoundary(Up)));
	procpoints.push_back(Vector2D(obc->GetGridBoundary(Left), obc->GetGridBoundary(Up)));

	vector<Vector2D> points = UpdatePoints(pv, obc);
	vector<int> HilbertIndeces;
	if (reorder)
	{
		HilbertIndeces = HilbertOrder(points, static_cast<int>(pv.size()));
		points = VectorValues(points, HilbertIndeces);
	}

	Tri.update(points, procpoints);

	Nextra = static_cast<int>(Tri.ChangeCor().size());
	vector<Edge> box_edges = obc->GetBoxEdges();
	Tri.BuildBoundary(obc, box_edges);

	eps = 1e-8;
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	build_v();

	if (logger)
		logger->output(*this);

	CM.resize(Tri.getCor().size());
	for (size_t i = 0; i < points.size(); ++i)
		CM[i] = CalcCellCM(i);

	size_t counter = pv.size() + 3;
	if (obc->GetBoundaryType() == Periodic)
	{
		for (int i = static_cast<int>(counter); i < Tri.GetCorSize(); ++i)
		{
			int NorgIndex = Tri.GetOrgIndex(i);
			if (NorgIndex < Nextra)
			{
			  CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + (Tri.get_point(static_cast<size_t>(i)) -
					Tri.get_point(static_cast<size_t>(NorgIndex)));
			}
		}
	}
	else
	{
		if (obc->GetBoundaryType() == Rectengular)
		{
			for (int i = static_cast<int>(counter); i < Tri.GetCorSize(); ++i)
			{
				int NorgIndex = Tri.GetOrgIndex(i);
				if (NorgIndex < Nextra)
				{
				  Vector2D norm = normalize(Tri.get_point(static_cast<size_t>(i)) - Tri.get_point(static_cast<size_t>(NorgIndex)));
				  CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + Tri.get_point(static_cast<size_t>(i))
						- Tri.get_point(static_cast<size_t>(NorgIndex)) - 2 * norm * ScalarProd(norm,
							CM[static_cast<size_t>(NorgIndex)] - Tri.get_point(static_cast<size_t>(NorgIndex)));
				}
			}
		}
		else  // Half periodic case
		{
			double dx = obc->GetGridBoundary(Right) - obc->GetGridBoundary(Left);
			for (int i = static_cast<int>(counter); i < Tri.GetCorSize(); ++i)
			{
				int NorgIndex = Tri.GetOrgIndex(i);
				if (NorgIndex < Nextra)
				{
				  double dx_temp = fastabs(Tri.get_point(static_cast<size_t>(i)) - Tri.get_point(static_cast<size_t>(NorgIndex)));
					if (dx_temp<1.0001*dx && dx_temp*1.0001>dx)
					  CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + (Tri.get_point(static_cast<size_t>(i)) -
							Tri.get_point(static_cast<size_t>(NorgIndex)));
					else
					{
					  Vector2D norm = normalize(Tri.get_point(static_cast<size_t>(i)) - Tri.get_point(static_cast<size_t>(NorgIndex)));
					  CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + Tri.get_point(static_cast<size_t>(i))
							- Tri.get_point(static_cast<size_t>(NorgIndex)) - 2 * norm * ScalarProd(norm,
								CM[static_cast<size_t>(NorgIndex)] - Tri.get_point(static_cast<size_t>(NorgIndex)));
					}
				}
			}
		}
	}
	return HilbertIndeces;
}

vector<int> VoronoiMesh::GetNeighbors(int index)const
{
	vector<int> res(mesh_vertices[static_cast<size_t>(index)].size());
	for (size_t i = 0; i < res.size(); ++i)
		res[i] = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][i])].neighbors.first != index ?
		edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][i])].neighbors.first :
		edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][i])].neighbors.second;
	return res;
}

void VoronoiMesh::GetNeighbors(int index, vector<int> &neigh)const
{
	neigh.resize(mesh_vertices[static_cast<size_t>(index)].size());
	size_t N = neigh.size();
	for (size_t i = 0; i < N; ++i)
		neigh[i] = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][i])].neighbors.first != index ?
		edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][i])].neighbors.first :
		edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][i])].neighbors.second;
}


vector<int> VoronoiMesh::GetLiteralNeighbors(int index)const
{
	int n = static_cast<int>(mesh_vertices[static_cast<size_t>(index)].size());
	vector<int> res;
	res.reserve(static_cast<size_t>(n));
	for (int i = 0; i < n; ++i)
	{
		int other = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].neighbors.first;
		if (other != index)
		{
			if (other > -1)
				res.push_back(other);
		}
		else
		{
			if (other > -1)
				other = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(index)][static_cast<size_t>(i)])].neighbors.second;
			res.push_back(other);
		}
	}
	return res;
}

VoronoiMesh* VoronoiMesh::clone(void)const
{
	return new VoronoiMesh(*this);
}

int VoronoiMesh::GetPointNo(void) const
{
	return Tri.get_length();
}

void VoronoiMesh::SetPointNo(int N)
{
	Tri.Changelength(N);
	return;
}

Vector2D VoronoiMesh::GetMeshPoint(int index) const
{
	return Tri.get_point(static_cast<size_t>(index));
}

int VoronoiMesh::GetTotalSidesNumber(void) const
{
	return static_cast<int>(edges.size());
}

const vector<Edge>& VoronoiMesh::getAllEdges(void) const
{
	return edges;
}

Edge const& VoronoiMesh::GetEdge(int index) const
{
	return edges[static_cast<size_t>(index)];
}

Vector2D const& VoronoiMesh::GetCellCM(int index) const
{
	return CM[static_cast<size_t>(index)];
}

vector<Vector2D>& VoronoiMesh::GetAllCM(void)
{
	return CM;
}

void VoronoiMesh::FindIntersectingOuterPoints(vector<Edge> const&box_edges, vector<vector<int> >
	&boxduplicate, vector<vector<int> > const&firstduplicated)
{
	int n = static_cast<int>(box_edges.size());
	boxduplicate.resize(static_cast<size_t>(n));
	int N = static_cast<int>(mesh_vertices.size());
	if (N < 20)
	{
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < N; ++j)
				boxduplicate[static_cast<size_t>(i)].push_back(j);
		return;
	}

	N = static_cast<int>(firstduplicated.size());
	for (int i = 0; i < N; ++i)
	{
		n = static_cast<int>(firstduplicated[static_cast<size_t>(i)].size());
		for (int j = 0; j < n; ++j)
		{
			vector<int> temp = CellIntersectOuterBoundary(box_edges, firstduplicated[static_cast<size_t>(i)][static_cast<size_t>(j)]);
			int jj = static_cast<int>(temp.size());
			if (jj > 0)
			{
				for (int k = 0; k < jj; ++k)
					boxduplicate[static_cast<size_t>(temp[static_cast<size_t>(k)])].push_back(firstduplicated[static_cast<size_t>(i)][static_cast<size_t>(j)]);
			}
		}
	}
	n = static_cast<int>(box_edges.size());
	for (int i = 0; i < n; ++i)
	{
		if (!boxduplicate[static_cast<size_t>(i)].empty())
		{
			sort(boxduplicate[static_cast<size_t>(i)].begin(), boxduplicate[static_cast<size_t>(i)].end());
			boxduplicate[static_cast<size_t>(i)] = unique(boxduplicate[static_cast<size_t>(i)]);
		}
	}
}

void VoronoiMesh::FindIntersectingPoints(vector<Edge> const &box_edges,
	vector<vector<int> > &toduplicate)
{
	int n = static_cast<int>(box_edges.size());
	toduplicate.resize(static_cast<size_t>(n));
	int N = static_cast<int>(mesh_vertices.size());
	if (N < 20)
	{
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < N; ++j)
				toduplicate[static_cast<size_t>(i)].push_back(j);
		return;
	}

	for (int i = 0; i < N; ++i)
	{
		vector<int> temp;
		temp = CellIntersectBoundary(box_edges, i);
		int j = static_cast<int>(temp.size());
		if (j > 0)
		{
			for (int k = 0; k < j; ++k)
				toduplicate[static_cast<size_t>(temp[static_cast<size_t>(k)])].push_back(i);
		}
	}
	for (int i = 0; i < n; ++i)
	{
		if (!toduplicate[static_cast<size_t>(i)].empty())
		{
			sort(toduplicate[static_cast<size_t>(i)].begin(), toduplicate[static_cast<size_t>(i)].end());
			toduplicate[static_cast<size_t>(i)] = unique(toduplicate[static_cast<size_t>(i)]);
		}
	}
}

vector<int> VoronoiMesh::CellIntersectBoundary(vector<Edge> const&box_edges, int cell)
{
	int ncell = static_cast<int>(mesh_vertices[static_cast<size_t>(cell)].size());
	int nbox = static_cast<int>(box_edges.size());
	vector<int> res;
	Vector2D intersect;
	for (int i = 0; i < ncell; ++i)
	{
		for (int j = 0; j < nbox; ++j)
		{
			if (SegmentIntersection(box_edges[static_cast<size_t>(j)], edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])],
				intersect))
				res.push_back(j);
		}
	}
	sort(res.begin(), res.end());
	res = unique(res);
	int nintersect = static_cast<int>(res.size());
	if (nintersect > 1)
	{
		vector<Vector2D> cpoints;
		ConvexHull(cpoints, *this, cell);
		for (int i = 0; i < nbox; ++i)
			if (PointInCell(cpoints, box_edges[static_cast<size_t>(i)].vertices.first) ||
				PointInCell(cpoints, box_edges[static_cast<size_t>(i)].vertices.second))
				res.push_back(i);
		sort(res.begin(), res.end());
		res = unique(res);
	}
	return res;
}

vector<int> VoronoiMesh::CellIntersectOuterBoundary(vector<Edge> const&box_edges, int cell)
{
	int ncell = static_cast<int>(mesh_vertices[static_cast<size_t>(cell)].size());
	int nbox = static_cast<int>(box_edges.size());
	vector<int> res;
	Vector2D intersect;
	boost::array<Vector2D, 3> tocheck;
	for (int i = 0; i < ncell; ++i)
	{
		for (int j = 0; j < nbox; ++j)
		{
			if (SegmentIntersection(box_edges[static_cast<size_t>(j)], edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])],
				intersect))
			{
				double r = sqrt(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].GetLength()*
					box_edges[static_cast<size_t>(j)].GetLength());
				double eps1 = 1e-7;
				if (abs(orient2d(TripleConstRef<Vector2D>
					(box_edges[static_cast<size_t>(j)].vertices.second - box_edges[static_cast<size_t>(j)].vertices.first,
						edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.second - box_edges[static_cast<size_t>(j)].vertices.first,
						edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.first - box_edges[static_cast<size_t>(j)].vertices.first)))
					< r*r*eps1)
					continue;
				if (DistanceToEdge(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.first,
					box_edges[static_cast<size_t>(j)]) < eps1*r)
					continue;
				if (DistanceToEdge(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.second,
					box_edges[static_cast<size_t>(j)]) < eps1*r)
					continue;
				if ((intersect.distance(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.first)
					> eps1*r) && (intersect.distance(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(cell)][static_cast<size_t>(i)])].vertices.second)
					> eps1*r))
					res.push_back(j);
			}
		}
	}
	sort(res.begin(), res.end());
	res = unique(res);
	return res;
}


bool VoronoiMesh::CloseToBorder(int point, int &border)
{
	int olength = Tri.GetOriginalLength();
	int n = static_cast<int>(mesh_vertices[static_cast<size_t>(point)].size());
	for (int i = 0; i < n; ++i)
	{
		if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second == point)
			border = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first;
		else
			border = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second;
		if (border > olength)
			return true;
	}
	return false;
}

vector<int> VoronoiMesh::GetBorderingCells(vector<int> const& copied,
	vector<int> const& totest, int tocheck, vector<int> tempresult, int outer)
{
	int border, test;
	int olength = Tri.GetOriginalLength();
	tempresult.push_back(tocheck);
	sort(tempresult.begin(), tempresult.end());
	int n = static_cast<int>(mesh_vertices[static_cast<size_t>(tocheck)].size());
	for (int i = 0; i < n; ++i)
	{
		if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(tocheck)][static_cast<size_t>(i)])].neighbors.second == tocheck)
			test = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(tocheck)][static_cast<size_t>(i)])].neighbors.first;
		else
			test = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(tocheck)][static_cast<size_t>(i)])].neighbors.second;
		if (test >= olength)
			continue;
		if (test < 0)
			continue;
		if (CloseToBorder(test, border))
			if (border == outer)
				if (!binary_search(copied.begin(), copied.end(), test) &&
					!binary_search(totest.begin(), totest.end(), test) &&
					!binary_search(tempresult.begin(), tempresult.end(), test))
					tempresult = GetBorderingCells(copied, totest, test, tempresult, outer);
	}
	return tempresult;
}

void VoronoiMesh::GetAdditionalBoundary(vector<vector<int> > &copied,
	vector<vector<int> > &neighbors, vector<vector<int> > &totest)
{
	int nsides = static_cast<int>(copied.size());
	// Get all the neighbors
	neighbors.clear();
	neighbors.resize(static_cast<size_t>(nsides));
	for (int i = 0; i < nsides; ++i)
	{
		sort(copied[static_cast<size_t>(i)].begin(), copied[static_cast<size_t>(i)].end());
		// look if there are boundary points neighbors
		int n = static_cast<int>(totest[static_cast<size_t>(i)].size());
		for (int j = 0; j < n; ++j)
		{
			if (totest[static_cast<size_t>(i)][static_cast<size_t>(j)] == -1)
				continue;
			vector<int> toadd;
			int outer = 0;
			if (CloseToBorder(totest[static_cast<size_t>(i)][static_cast<size_t>(j)], outer))
				toadd = GetBorderingCells(copied[static_cast<size_t>(i)], totest[static_cast<size_t>(i)], totest[static_cast<size_t>(i)][static_cast<size_t>(j)], toadd, outer);
			int nn = static_cast<int>(toadd.size());
			for (int k = 0; k < nn; ++k)
				neighbors[static_cast<size_t>(i)].push_back(toadd[static_cast<size_t>(k)]);
		}
		sort(neighbors[static_cast<size_t>(i)].begin(), neighbors[static_cast<size_t>(i)].end());
		neighbors[static_cast<size_t>(i)] = unique(neighbors[static_cast<size_t>(i)]);
		neighbors[static_cast<size_t>(i)] = RemoveList(neighbors[static_cast<size_t>(i)], copied[static_cast<size_t>(i)]);
	}
}

void VoronoiMesh::GetRealNeighbor(vector<int> &result, int point) const
{
	result.reserve(7);
	int n = static_cast<int>(mesh_vertices[static_cast<size_t>(point)].size());
	int olength = Tri.GetOriginalLength();
	for (int i = 0; i < n; ++i)
	{
		if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first == point)
		{
			if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second > -1 &&
				edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second < olength)
				result.push_back(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.second);
		}
		else
		{
			if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first > -1 &&
				edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first < olength)
				result.push_back(edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].neighbors.first);
		}
	}
	sort(result.begin(), result.end());
	result = unique(result);
}

void VoronoiMesh::GetNeighborNeighbors(vector<int> &result, int point) const
{
	vector<int> neigh = GetNeighbors(point);
	result.clear();
	result.reserve(25);
	//GetRealNeighbor(neigh,point);
	int n = static_cast<int>(neigh.size());
	for (int i = 0; i < n; ++i)
	{
#ifdef RICH_MPI
		if (GetOriginalIndex(neigh[static_cast<size_t>(i)]) > GetPointNo())
			continue;
#endif
		vector<int> temp = GetNeighbors(GetOriginalIndex(neigh[static_cast<size_t>(i)]));
		for (size_t j = 0; j < temp.size(); ++j)
			result.push_back(GetOriginalIndex(temp[j]));
	}
	sort(result.begin(), result.end());
	result = unique(result);
	// Remove self point
	RemoveVal(result, point);
}

void VoronoiMesh::GetNeighborNeighborsMPI(vector<int> &result, int point)
{
	vector<int> neigh;
	result.clear();
	GetRealNeighbor(neigh, point);
	int n = static_cast<int>(neigh.size());
	for (int i = 0; i < n; ++i)
	{
		vector<int> temp;
		GetRealNeighbor(temp, neigh[static_cast<size_t>(i)]);
		int N = static_cast<int>(temp.size());
		for (int j = 0; j < N; ++j)
			result.push_back(temp[static_cast<size_t>(j)]);
	}
	sort(result.begin(), result.end());
	unique(result);
}

void VoronoiMesh::GetCorners(vector<vector<int> > &copied,
	vector<vector<int> > &result)
{
	// copied should be sorted already
	int nsides = static_cast<int>(copied.size());
	result.clear();
	OrgCorner.clear();
	OrgCorner.resize(static_cast<size_t>(nsides));
	result.resize(static_cast<size_t>(nsides));
	vector<vector<int> > toadd(static_cast<size_t>(nsides));
	for (int i = 0; i < nsides; ++i)
	{
		int n = static_cast<int>(copied[static_cast<size_t>(i)].size());
		for (int j = 0; j < n; ++j)
		{
			if (binary_search(copied[static_cast<size_t>((i + 1) % nsides)].begin(), copied[static_cast<size_t>((i + 1) % nsides)].end(),
				copied[static_cast<size_t>(i)][static_cast<size_t>(j)]))
			{
				vector<int> temp;
				GetNeighborNeighborsMPI(temp, copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
				result[static_cast<size_t>(i)].insert(result[static_cast<size_t>(i)].end(), temp.begin(), temp.end());
				temp = AddPointsAlongEdge(static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]), copied, i);
				toadd[static_cast<size_t>((i + 1) % nsides)].insert(toadd[static_cast<size_t>((i + 1) % nsides)].end(), temp.begin(),
					temp.end());
				temp = GetNeighbors(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
				for (vector<int>::iterator it = temp.begin(); it != temp.end(); ++it)
					if (*it<GetPointNo() && *it>-1)
						OrgCorner[static_cast<size_t>(i)].push_back(*it);
				OrgCorner[static_cast<size_t>(i)].push_back(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
			}
			if (binary_search(copied[static_cast<size_t>((i - 1 + nsides) % nsides)].begin(), copied[static_cast<size_t>((i - 1 + nsides) % nsides)].end(),
				copied[static_cast<size_t>(i)][static_cast<size_t>(j)]))
			{
				vector<int> temp;
				GetNeighborNeighborsMPI(temp, copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
				result[static_cast<size_t>((i - 1 + nsides) % nsides)].insert(result[static_cast<size_t>((i - 1 + nsides) % nsides)].end()
					, temp.begin(), temp.end());
				temp = AddPointsAlongEdge(static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]), copied, i);
				toadd[static_cast<size_t>((i - 1 + nsides) % nsides)].insert(toadd[static_cast<size_t>((i - 1 + nsides) % nsides)].end(),
					temp.begin(), temp.end());
				temp = GetNeighbors(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
				for (vector<int>::iterator it = temp.begin(); it != temp.end(); ++it)
					if (*it<GetPointNo() && *it>-1)
						OrgCorner[static_cast<size_t>((i - 1 + nsides) % nsides)].push_back(*it);
				OrgCorner[static_cast<size_t>((i - 1 + nsides) % nsides)].push_back(copied[static_cast<size_t>(i)][static_cast<size_t>(j)]);
			}
		}
	}
	for (int i = 0; i < nsides; ++i)
	{
		copied[static_cast<size_t>(i)].insert(copied[static_cast<size_t>(i)].end(), toadd[static_cast<size_t>(i)].begin(), toadd[static_cast<size_t>(i)].end());
		sort(copied[static_cast<size_t>(i)].begin(), copied[static_cast<size_t>(i)].end());
		copied[static_cast<size_t>(i)] = unique(copied[static_cast<size_t>(i)]);
		sort(result[static_cast<size_t>(i)].begin(), result[static_cast<size_t>(i)].end());
		result[static_cast<size_t>(i)] = unique(result[static_cast<size_t>(i)]);
		if (!OrgCorner[static_cast<size_t>(i)].empty())
		{
			sort(OrgCorner[static_cast<size_t>(i)].begin(), OrgCorner[static_cast<size_t>(i)].end());
			OrgCorner[static_cast<size_t>(i)] = unique(OrgCorner[static_cast<size_t>(i)]);
		}
	}
}

void VoronoiMesh::GetToTest(vector<vector<int> > &copied, vector<vector<int> > &totest)
{
	int nsides = static_cast<int>(copied.size());
	int olength = Tri.GetOriginalLength();
	// sort the vectors
	for (int i = 0; i < nsides; ++i)
		sort(copied[static_cast<size_t>(i)].begin(), copied[static_cast<size_t>(i)].end());
	totest.resize(static_cast<size_t>(nsides));
	int test = 0;
	for (int i = 0; i < nsides; ++i)
	{
		vector<int> totest2;
		int ncopy = static_cast<int>(copied[static_cast<size_t>(i)].size());
		for (int j = 0; j < ncopy; ++j)
		{
			int n = static_cast<int>(mesh_vertices[static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)])].size());
			for (int k = 0; k < n; ++k)
			{
				if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)])][static_cast<size_t>(k)])].neighbors.first ==
					copied[static_cast<size_t>(i)][static_cast<size_t>(j)])
					test = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)])][static_cast<size_t>(k)])].neighbors.second;
				else
					test = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(copied[static_cast<size_t>(i)][static_cast<size_t>(j)])][static_cast<size_t>(k)])].neighbors.first;
				if (test < olength)
					totest2.push_back(test);
			}
		}
		sort(totest2.begin(), totest2.end());
		totest2 = unique(totest2);
		totest[static_cast<size_t>(i)] = totest2;
	}
}

vector<int> VoronoiMesh::FindEdgeStartConvex(int point)
{
	int n = static_cast<int>(mesh_vertices[static_cast<size_t>(point)].size());
	Vector2D min_point;
	int min_index = 0, p_index;
	if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(0)])].vertices.first.x <
		edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][0])].vertices.second.x)
	{
		min_point = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][0])].vertices.first;
		p_index = 0;
	}
	else
	{
		min_point = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][0])].vertices.second;
		p_index = 1;
	}
	for (int i = 1; i < n; ++i)
	{
		double R = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].GetLength();
		if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first.x < (min_point.x - R*eps))
		{
			min_point = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first;
			min_index = i;
			p_index = 0;
		}
		if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second.x < (min_point.x - R*eps))
		{
			min_point = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second;
			min_index = i;
			p_index = 1;
		}
		if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first.x < (min_point.x + R*eps) &&
			edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first.y < min_point.y)
		{
			min_point = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.first;
			min_index = i;
			p_index = 0;
		}

		if (edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second.x < (min_point.x + R*eps) &&
			edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second.y < min_point.y)
		{
			min_point = edges[static_cast<size_t>(mesh_vertices[static_cast<size_t>(point)][static_cast<size_t>(i)])].vertices.second;
			min_index = i;
			p_index = 1;
		}
	}
	vector<int> res(2);
	res[0] = min_index;
	res[1] = p_index;
	return res;
}

void VoronoiMesh::ConvexEdgeOrder(void)
{
	int n = static_cast<int>(mesh_vertices.size());
	for (int i = 0; i < n; ++i)
	{
		double R = GetWidth(i);
		vector<int> min_index = FindEdgeStartConvex(i);
		int p_loc = min_index[1];
		int edge_loc = mesh_vertices[static_cast<size_t>(i)][static_cast<size_t>(min_index[0])];
		int nedges = static_cast<int>(mesh_vertices[static_cast<size_t>(i)].size());
		std::list<int> elist;
		for (int j = 0; j < nedges; ++j)
		{
			if (j != min_index[0])
				elist.push_back(mesh_vertices[static_cast<size_t>(i)][static_cast<size_t>(j)]);
		}
		vector<int> new_order;
		new_order.reserve(static_cast<size_t>(nedges));
		new_order.push_back(edge_loc);
		for (int j = 0; j < nedges; ++j)
		{
			if (j == min_index[0])
				continue;
			int nlist = static_cast<int>(elist.size());
			std::list<int>::iterator it = elist.begin();
			for (int k = 0; k < nlist; ++k)
			{
				double temp0 = pair_member(edges[static_cast<size_t>(edge_loc)].vertices, (p_loc + 1) % 2).distance(edges[static_cast<size_t>(*it)].vertices.first);
				if (temp0 < eps*R)
				{
					p_loc = 0;
					edge_loc = *it;
					elist.erase(it);
					new_order.push_back(edge_loc);
					break;
				}
				double temp1 = pair_member(edges[static_cast<size_t>(edge_loc)].vertices, (p_loc + 1) % 2).distance(
					edges[static_cast<size_t>(*it)].vertices.second);
				if (temp1 < eps*R)
				{
					p_loc = 1;
					edge_loc = *it;
					new_order.push_back(edge_loc);
					elist.erase(it);
					break;
				}
				++it;
			}
		}
		mesh_vertices[static_cast<size_t>(i)] = new_order;
	}
}

vector<Edge>& VoronoiMesh::GetAllEdges(void)
{
	return edges;
}

vector<vector<int> >& VoronoiMesh::GetDuplicatedPoints(void)
{
	return GhostPoints;
}

vector<vector<int> >const& VoronoiMesh::GetDuplicatedPoints(void)const
{
	return GhostPoints;
}

vector<vector<int> >& VoronoiMesh::GetGhostIndeces(void)
{
	return NGhostReceived;
}

vector<vector<int> >const& VoronoiMesh::GetGhostIndeces(void)const
{
	return NGhostReceived;
}

int VoronoiMesh::GetTotalPointNumber(void)const
{
	return static_cast<int>(Tri.getCor().size());
}

vector<int> VoronoiMesh::GetDuplicatedProcs(void)const
{
	return GhostProcs;
}

vector<int> VoronoiMesh::GetSentProcs(void)const
{
	return SentProcs;
}

vector<vector<int> >const& VoronoiMesh::GetSentPoints(void)const
{
	return SentPoints;
}

// cpoints must be convex hull, checks if vec is inside cpoints
bool PointInCell(vector<Vector2D> const& cpoints, Vector2D const& vec)
{
	for (size_t i = 0, endp = cpoints.size(); i < endp; ++i)
	{
		if (orient2d(TripleConstRef<Vector2D>(cpoints[i],
			cpoints[(i + 1) % endp],
			vec)) < 0)
			return false;
	}
	return true;
}

// result is : minx, maxx, miny, maxy
boost::array<double, 4> VoronoiMesh::FindMaxCellEdges(void)
{
	int n = static_cast<int>(cell_edges.size());
	boost::array<double, 4> res;
	res[0] = min(cell_edges[0].vertices.first.x, cell_edges[0].vertices.second.x);
	res[1] = max(cell_edges[0].vertices.first.x, cell_edges[0].vertices.second.x);
	res[2] = min(cell_edges[0].vertices.first.y, cell_edges[0].vertices.second.y);
	res[3] = max(cell_edges[0].vertices.first.y, cell_edges[0].vertices.second.y);
	for (int i = 1; i < n; ++i)
	{
		res[0] = min(min(cell_edges[static_cast<size_t>(i)].vertices.first.x, cell_edges[static_cast<size_t>(i)].vertices.second.x), res[0]);
		res[1] = max(max(cell_edges[static_cast<size_t>(i)].vertices.first.x, cell_edges[static_cast<size_t>(i)].vertices.second.x), res[1]);
		res[2] = min(min(cell_edges[static_cast<size_t>(i)].vertices.first.y, cell_edges[static_cast<size_t>(i)].vertices.second.y), res[2]);
		res[3] = max(max(cell_edges[static_cast<size_t>(i)].vertices.first.y, cell_edges[static_cast<size_t>(i)].vertices.second.y), res[3]);
	}
	return res;
}

#ifdef RICH_MPI
vector<int> VoronoiMesh::Update
(const vector<Vector2D>& points,
	const Tessellation& vproc, bool reorder)
{
	NGhostReceived.clear();
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<int> cedges;
	ConvexEdges(cedges, vproc, rank);
	cell_edges.clear();
	for (size_t i = 0; i < cedges.size(); ++i)
		cell_edges.push_back(vproc.GetEdge(cedges[i]));

	// Get the convex hull of the cell
	vector<Vector2D> cpoints;
	ConvexHull(cpoints, vproc, rank);

	vector<int> HilbertIndeces;
	vector<Vector2D> newcor;
	if (reorder)
	{
		HilbertIndeces = HilbertOrder(points, static_cast<int>(points.size()));
		newcor = VectorValues(points, HilbertIndeces);
	}
	else
		newcor = points;


	// Did points move between procs?
	newcor = UpdateMPIPoints(vproc, rank, newcor, obc, selfindex, SentProcs, SentPoints);

	//Build the delaunay
	try
	{
		Tri.build_delaunay(newcor, cpoints);
	}
	catch (UniversalError & eo)
	{
		size_t bad_index = static_cast<size_t>(eo.GetValues()[0]);
		eo.AddEntry("rank", rank);
		eo.AddEntry("original point in cor x", points[bad_index].x);
		eo.AddEntry("original point in cor y", points[bad_index].y);
		vector<int> edge_index = vproc.GetCellEdges(rank);
		for (size_t i = 0; i < edge_index.size(); ++i)
		{
			eo.AddEntry("edge number", static_cast<double>(edge_index[i]));
			eo.AddEntry("Edge v0x", vproc.GetEdge(edge_index[i]).vertices.first.x);
			eo.AddEntry("Edge v0y", vproc.GetEdge(edge_index[i]).vertices.first.y);
			eo.AddEntry("Edge v1x", vproc.GetEdge(edge_index[i]).vertices.second.x);
			eo.AddEntry("Edge v1y", vproc.GetEdge(edge_index[i]).vertices.second.y);
		}
		throw eo;
	}
	Nextra = static_cast<int>(Tri.ChangeCor().size());
	eps = 1e-8;
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	NGhostReceived.clear();
	pair<vector<vector<int> >, vector<int> > ptemp = Tri.BuildBoundary(obc, vproc, NGhostReceived);
	GhostPoints = ptemp.first;
	GhostProcs = ptemp.second;
	build_v();

	if (logger)
		logger->output(*this);

	size_t n = static_cast<size_t>(GetPointNo());
	CM.resize(Tri.getCor().size());
	for (size_t i = 0; i < n; ++i)
		CM[i] = CalcCellCM(i);

	if (obc->GetBoundaryType() == Periodic)
	{
		for (int i = Nextra; i < Tri.GetCorSize(); ++i)
		{
			int NorgIndex = Tri.GetOrgIndex(i);
			if (NorgIndex < Nextra)
			{
				CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + (Tri.get_point(i) -
					Tri.get_point(static_cast<size_t>(NorgIndex)));
			}
		}
	}
	else
	{
		if (obc->GetBoundaryType() == Rectengular)
		{
			for (int i = Nextra; i < Tri.GetCorSize(); ++i)
			{
				int NorgIndex = Tri.GetOrgIndex(i);
				if (NorgIndex < Nextra)
				{
					Vector2D norm = normalize(Tri.get_point(i) - Tri.get_point(static_cast<size_t>(NorgIndex)));
					CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + Tri.get_point(i)
						- Tri.get_point(static_cast<size_t>(NorgIndex)) - 2 * norm * ScalarProd(norm,
							CM[static_cast<size_t>(NorgIndex)] - Tri.get_point(static_cast<size_t>(NorgIndex)));
				}
			}
		}
		else  // Half periodic case
			throw(UniversalError("No HalfPeriodic in MPI"));
	}
	// communicate the ghost CM
	vector<vector<Vector2D> > incoming = MPI_exchange_data(GhostProcs, GhostPoints, CM);
	// Add the recieved CM
	for (size_t i = 0; i < incoming.size(); ++i)
		for (size_t j = 0; j < incoming.at(i).size(); ++j)
			CM[NGhostReceived.at(i).at(j)] = incoming[i][j];
	// Fix for periodic changes
	if (obc->GetBoundaryType() == Periodic)
	{
		std::vector<Vector2D> & localcor = Tri.ChangeCor();
		incoming = MPI_exchange_data(GhostProcs, GhostPoints, localcor);
		for (size_t i = 0; i < incoming.size(); ++i)
			for (size_t j = 0; j < incoming.at(i).size(); ++j)
				CM[NGhostReceived.at(i).at(j)] += localcor[NGhostReceived.at(i).at(j)] - incoming[i][j];
	}

	return HilbertIndeces;
}

void VoronoiMesh::Initialise
(vector<Vector2D> const& pv,
	Tessellation const& vproc,
	OuterBoundary const* outer, bool reorder)
{
	NGhostReceived.clear();
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	obc = outer;
	vector<int> cedges;
	ConvexEdges(cedges, vproc, rank);
	cell_edges.clear();
	for (size_t i = 0; i < cedges.size(); ++i)
		cell_edges.push_back(vproc.GetEdge(cedges[i]));

	vector<Vector2D> points;
	if (reorder)
		points = VectorValues(pv, HilbertOrder(pv, static_cast<int>(pv.size())));
	else
		points = pv;

	// Get the convex hull of the cell
	vector<Vector2D> cpoints;
	ConvexHull(cpoints, vproc, rank);
	//Build the delaunay
	Tri.build_delaunay(points, cpoints);
	Nextra = static_cast<int>(Tri.ChangeCor().size());
	eps = 1e-8;
	edges.clear();
	GhostPoints.clear();
	GhostProcs.clear();
	NGhostReceived.clear();
	selfindex.resize(points.size());
	size_t npoints = points.size();
	for (size_t i = 0; i < npoints; ++i)
		selfindex[i] = i;
	pair<vector<vector<int> >, vector<int> > ptemp = Tri.BuildBoundary(outer, vproc, NGhostReceived);
	GhostPoints = ptemp.first;
	GhostProcs = ptemp.second;

	/*
	if(get_rank()==0){
	  {
	    ofstream dump("GhostProcs.txt");
	    for(size_t i=0;i<GhostProcs.size();++i)
	      dump << GhostProcs.at(i) << endl;
	  dump.close();
	  }
	  {
	    ofstream dump("GhostPoints.txt");
	    for(size_t i=0;i<GhostPoints.size();++i){
	      for(size_t j=0;j<GhostPoints.at(i).size();++j)
		dump << GhostPoints.at(i).at(j) << " ";
	      dump << endl;
	    }
	    dump.close();
	  }
	  assert(false);
	}
	*/

	build_v();

	if (logger)
		logger->output(*this);

	size_t n = static_cast<size_t>(GetPointNo());
	CM.resize(Tri.getCor().size());
	for (size_t i = 0; i < n; ++i)
		CM[i] = CalcCellCM(i);

	if (obc->GetBoundaryType() == Periodic)
	{
		for (int i = Nextra; i < Tri.GetCorSize(); ++i)
		{
			int NorgIndex = Tri.GetOrgIndex(i);
			if (NorgIndex < Nextra)
			{
				CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + (Tri.get_point(i) -
					Tri.get_point(static_cast<size_t>(NorgIndex)));
			}
		}
	}
	else
	{
		if (obc->GetBoundaryType() == Rectengular)
		{
			for (int i = Nextra; i < Tri.GetCorSize(); ++i)
			{
				int NorgIndex = Tri.GetOrgIndex(i);
				if (NorgIndex < Nextra)
				{
					Vector2D norm = normalize(Tri.get_point(i) - Tri.get_point(static_cast<size_t>(NorgIndex)));
					CM[static_cast<size_t>(i)] = CM[static_cast<size_t>(NorgIndex)] + Tri.get_point(i)
						- Tri.get_point(static_cast<size_t>(NorgIndex)) - 2 * norm * ScalarProd(norm,
							CM[static_cast<size_t>(NorgIndex)] - Tri.get_point(static_cast<size_t>(NorgIndex)));
				}
			}
		}
		else  // Half periodic case
			throw(UniversalError("No HalfPeriodic in MPI"));
	}
	// communicate the ghost CM
	vector<vector<Vector2D> > incoming = MPI_exchange_data(GhostProcs, GhostPoints, CM);
	// Add the recieved CM
	for (size_t i = 0; i < incoming.size(); ++i)
		for (size_t j = 0; j < incoming.at(i).size(); ++j)
			CM[NGhostReceived.at(i).at(j)] = incoming[i][j];

	// Fix for periodic changes
	if (obc->GetBoundaryType() == Periodic)
	{
		std::vector<Vector2D> & localcor = Tri.ChangeCor();
		incoming = MPI_exchange_data(GhostProcs, GhostPoints, localcor);
		for (size_t i = 0; i < incoming.size(); ++i)
			for (size_t j = 0; j < incoming.at(i).size(); ++j)
				CM[NGhostReceived.at(i).at(j)] += localcor[NGhostReceived.at(i).at(j)] - incoming[i][j];
	}
}


vector<Vector2D> VoronoiMesh::UpdateMPIPoints(Tessellation const& vproc, int rank,
	vector<Vector2D> &points, OuterBoundary const* obc, vector<size_t> &selfindex,
	vector<int> &sentproc, vector<vector<int> > &sentpoints)
{
	vector<Vector2D> res;
	res.reserve(points.size());
	selfindex.clear();
	size_t npoints = points.size();
	size_t nproc = static_cast<size_t>(vproc.GetPointNo());
	const double dx = obc->GetGridBoundary(Right) - obc->GetGridBoundary(Left);
	const double dy = obc->GetGridBoundary(Up) - obc->GetGridBoundary(Down);
	vector<Vector2D> cproc;
	ConvexHull(cproc, vproc, rank);
	vector<int> neighbors = vproc.GetNeighbors(rank);
	vector<int> realneigh;
	std::vector<size_t> neigh_keep;
	vector<vector<Vector2D> > neigh_chull;
	sentpoints.clear();
	sentproc.clear();
	bool periodic = obc->GetBoundaryType() == Periodic;
	for (size_t i = 0; i < neighbors.size(); ++i)
	{
		if (static_cast<size_t>(neighbors[i]) < nproc || periodic)
		{
			if (static_cast<size_t>(neighbors[i]) >= nproc)
			{
				int temp_neigh = vproc.GetOriginalIndex(neighbors[i]);
				if (temp_neigh == rank)
					continue;
				realneigh.push_back(temp_neigh);
				vector<Vector2D> temp;
				ConvexHull(temp, vproc, temp_neigh);
				Vector2D to_add_neigh = vproc.GetMeshPoint(neighbors[i]) - vproc.GetMeshPoint(temp_neigh);
				for (size_t j = 0; j < temp.size(); ++j)
					temp[j] += to_add_neigh;
				neigh_chull.push_back(temp);
				sentproc.push_back(temp_neigh);
				neigh_keep.push_back(i);
			}
			else
			{
				realneigh.push_back(neighbors[i]);
				vector<Vector2D> temp;
				ConvexHull(temp, vproc, neighbors[i]);
				neigh_chull.push_back(temp);
				sentproc.push_back(neighbors[i]);
				neigh_keep.push_back(i);
			}
		}
	}
	neighbors = VectorValues(neighbors, neigh_keep);
	std::sort(sentproc.begin(), sentproc.end());
	sentproc = unique(sentproc);
	sentpoints.resize(sentproc.size());

	for (size_t i = 0; i < npoints; ++i)
	{
		bool good = false;
		Vector2D temp = points[i];
		if (PointInCell(cproc, temp)) // Check own cpu
		{
			res.push_back(temp);
			selfindex.push_back(i);
			good = true;
		}
		if (!good)
		{
			for (size_t j = 0; j < realneigh.size(); ++j) // check cpu neighbors
			{
				if((!periodic && PointInCell(neigh_chull[j], temp)) || (periodic && PointInCell(neigh_chull[j], temp)))
				{
					good = true;
					if (periodic && static_cast<size_t>(neighbors[j]) >= nproc) // Do we need to move point?
						points[i] += vproc.GetMeshPoint(realneigh[j]) - vproc.GetMeshPoint(neighbors[j]);
					size_t index = std::find(sentproc.begin(), sentproc.end(), realneigh[j]) - sentproc.begin();
					assert(index < sentproc.size());
					sentpoints[index].push_back(static_cast<int>(i));
					break;
				}
			}
		}
		if (!good)
		{
			vector<Vector2D> cellpoints;
			for (size_t j = 0; j < nproc; ++j) // Search all cpus
			{
				if (std::find(realneigh.begin(), realneigh.end(), j) != realneigh.end() || j == static_cast<size_t>(rank))
					continue;
				ConvexHull(cellpoints, vproc, static_cast<int>(j));
				if ((!periodic && PointInCell(cellpoints, temp)) || (periodic && PointInCell(cellpoints, temp)))
				{
					good = true;
					size_t index = std::find(sentproc.begin(), sentproc.end(), j) - sentproc.begin();
					if (index >= sentproc.size())
					{
						sentproc.push_back(static_cast<int>(j));
						sentpoints.push_back(vector<int>(1, static_cast<int>(i)));
					}
					else
						sentpoints[index].push_back(static_cast<int>(i));
					break;
				}
			}
		}
		// If periodic check for all instances of duplication
		if (!good && obc->GetBoundaryType() == Periodic)
		{
			vector<Vector2D> cellpoints;
			for (size_t j = 0; j < nproc; ++j)
			{
				ConvexHull(cellpoints, vproc, static_cast<int>(j));
				// Create periodic instances of cpu
				std::vector<Vector2D> moved_point(8, temp);
				moved_point[0] += Vector2D(dx, 0);
				moved_point[1] += Vector2D(-dx, 0);
				moved_point[2] += Vector2D(dx, dy);
				moved_point[3] += Vector2D(dx, -dy);
				moved_point[4] += Vector2D(-dx, dy);
				moved_point[5] += Vector2D(-dx, -dy);
				moved_point[6] += Vector2D(0, dy);
				moved_point[7] += Vector2D(0, -dy);
				// Check if inside
				for (size_t k = 0; k < 8; ++k)
				{
					if (PointInCell(cellpoints, moved_point[k]))
					{
						good = true;
						points[i] = moved_point[k];
						if (static_cast<int>(j) == rank)
						{
							res.push_back(moved_point[k]);
							selfindex.push_back(i);
							break;
						}
						size_t index = std::find(sentproc.begin(), sentproc.end(), j) - sentproc.begin();
						if (index >= sentproc.size())
						{
							sentproc.push_back(static_cast<int>(j));
							sentpoints.push_back(vector<int>(1, static_cast<int>(i)));
						}
						else
							sentpoints[index].push_back(static_cast<int>(i));
						break;
					}
				}
			}
		}
		if (good)
			continue;
		HDF5Logger log("verror" + int2str(rank) + ".h5");
		log.output(vproc);
		UniversalError eo("Point is not inside any processor");
		eo.AddEntry("CPU rank", rank);
		eo.AddEntry("Point number", static_cast<double>(i));
		eo.AddEntry("Point x cor", points[i].x);
		eo.AddEntry("Point y cor", points[i].y);
		throw eo;
	}
	// Send/Recv the points
	// Communication
	int wsize;
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	vector<int> totalk(static_cast<size_t>(wsize), 0);
	vector<int> scounts(totalk.size(), 1);
	for (size_t i = 0; i < sentproc.size(); ++i)
		totalk[sentproc[i]] = 1;
	int nrecv;
	MPI_Reduce_scatter(&totalk[0], &nrecv, &scounts[0], MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);

	vector<MPI_Request> req(sentproc.size());
	for (size_t i = 0; i < sentproc.size(); ++i)
		MPI_Isend(&wsize, 1, MPI_INT, sentproc[i], 3, MPI_COMM_WORLD, &req[i]);
	vector<int> talkwithme;
	for (int i = 0; i < nrecv; ++i)
	{
		MPI_Status status;
		MPI_Recv(&wsize, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
		talkwithme.push_back(status.MPI_SOURCE);
	}
	if (req.size() > 0)
		MPI_Waitall(static_cast<int>(req.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	for (size_t i = 0; i < talkwithme.size(); ++i)
	{
		if (std::find(sentproc.begin(), sentproc.end(), talkwithme[i]) == sentproc.end())
		{
			sentproc.push_back(talkwithme[i]);
			sentpoints.push_back(vector<int>());
		}
	}
	// Point exchange
	vector<vector<Vector2D> > incoming(sentproc.size());
	vector<vector<double> > tosend(sentproc.size());
	vector<double> torecv;
	double dtemp = 0;
	Vector2D vtemp;
	req.clear();
	req.resize(sentproc.size());
	vector<int> indeces;
	vector<Vector2D> cortemphilbert;
	for (size_t i = 0; i < sentproc.size(); ++i)
	{
		const int dest = sentproc.at(i);
		if (!sentpoints.at(i).empty())
		{
			cortemphilbert = VectorValues(points, sentpoints.at(i));
			indeces = HilbertOrder(cortemphilbert, static_cast<int>(cortemphilbert.size()));
			tosend[i] = list_serialize(VectorValues(cortemphilbert, indeces));
			sentpoints[i] = VectorValues(sentpoints[i], indeces);
		}
		if (tosend[i].empty())
			MPI_Isend(&dtemp, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD, &req[i]);
		else
			MPI_Isend(&tosend[i][0], static_cast<int>(tosend[i].size()), MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &req[i]);
	}
	for (size_t i = 0; i < sentproc.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		torecv.resize(static_cast<size_t>(count));
		MPI_Recv(&torecv[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (status.MPI_TAG == 0)
		{
			size_t location = static_cast<size_t>(std::find(sentproc.begin(), sentproc.end(), status.MPI_SOURCE) -
				sentproc.begin());
			if (location >= sentproc.size())
				throw UniversalError("Bad location in mpi exchange");
			incoming[location] = list_unserialize(torecv, vtemp);
		}
		else
		{
			if (status.MPI_TAG != 1)
				throw UniversalError("Recv bad mpi tag");
		}
	}
	if (req.size() > 0)
		MPI_Waitall(static_cast<int>(req.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);

	// Combine the vectors
	for (size_t i = 0; i < incoming.size(); ++i)
		for (size_t j = 0; j < incoming[i].size(); ++j)
			res.push_back(incoming[i][j]);
	return res;
}



#endif 
