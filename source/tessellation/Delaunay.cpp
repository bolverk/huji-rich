#include "Delaunay.hpp"
#include <vector>
#include <cmath>
#include "../misc/triplet.hpp"
#include <boost/foreach.hpp>
#ifdef RICH_MPI
#include <mpi.h>
#endif // RICH_MPI

namespace {
	pair<int, int> find_diff(const facet& f1, const facet& f2)
	{
		if (f1.vertices.first != f2.vertices.first &&
			f1.vertices.first != f2.vertices.second &&
			f1.vertices.first != f2.vertices.third)
			return pair<int, int>(f1.vertices.first, 0);
		else if (f1.vertices.second != f2.vertices.first &&
			f1.vertices.second != f2.vertices.second &&
			f1.vertices.second != f2.vertices.third)
			return pair<int, int>(f1.vertices.second, 1);
		else if (f1.vertices.third != f2.vertices.first &&
			f1.vertices.third != f2.vertices.second &&
			f1.vertices.third != f2.vertices.third)
			return pair<int, int>(f1.vertices.third, 2);
		else
			throw UniversalError("Delaunay, Couldn't find difference bewteen two facets");
	}
}

Delaunay::DataOnlyForBuild::DataOnlyForBuild() :copied(vector<vector<char> >())
{}

Delaunay::DataOnlyForBuild::DataOnlyForBuild(DataOnlyForBuild const& other) :copied(other.copied) {}

Delaunay::DataOnlyForBuild& Delaunay::DataOnlyForBuild::operator=
(DataOnlyForBuild const& other)
{
	if (this != &other)
	{
		copied = other.copied;
	}
	return *this;
}

Delaunay::Delaunay(void) :
	lastFacet(0), CalcRadius(false),
	radius(vector<double>()), cell_points(vector<Vector2D>()),
	PointWasAdded(false),
	last_facet_added(0),
	f(vector<facet>()),
	cor(vector<Vector2D>()),
	length(0),
	olength(0), location_pointer(0), last_loc(0),
	logger(0)
#ifdef RICH_MPI
	,OrgIndex(vector<int>())
#endif
{}

Delaunay::Delaunay(Delaunay const& other) :
	lastFacet(other.lastFacet),
	CalcRadius(other.CalcRadius),
	radius(other.radius), cell_points(other.cell_points),
	PointWasAdded(other.PointWasAdded),
	last_facet_added(other.last_facet_added),
	f(other.f),
	cor(other.cor),
	length(other.length),
	olength(other.olength),
	location_pointer(other.location_pointer),
	last_loc(other.last_loc),
	logger(other.logger)
#ifdef RICH_MPI
	,OrgIndex(other.OrgIndex)
#endif
{}

Delaunay::~Delaunay(void)
{
	cor.clear();
	f.clear();
	cell_points.clear();
}

namespace
{
	// Checks if a point is inside a triangle
	bool InTriangle(const TripleConstRef<Vector2D>& tri,
		const Vector2D& point)
	{
		return (orient2d(TripleConstRef<Vector2D>(tri.first,
			tri.second,
			point)) > 0) &&
			(orient2d(TripleConstRef<Vector2D>(tri.second,
				tri.third,
				point)) > 0) &&
			(orient2d(TripleConstRef<Vector2D>(tri.third,
				tri.first,
				point)) > 0);
	}

	// Assume cell is orederd in convexhull counterclockwise
	bool InCell(vector<Vector2D> const& points, Vector2D const& p)
	{
		int n = static_cast<int>(points.size());
		for (int i = 0; i < n; ++i)
		{
			if (CrossProduct(points[static_cast<size_t>(i)] - p, points[static_cast<size_t>((i + 1) % n)] - p) < 0)
				return false;
		}
		return true;
	}

	vector<double> CellSize(vector<Vector2D> const& points)
	{
		int n = static_cast<int>(points.size());
		double minx = points[0].x;
		double miny = points[0].y;
		double maxx = minx;
		double maxy = miny;
		for (int i = 1; i < n; ++i)
		{
			minx = min(points[static_cast<size_t>(i)].x, minx);
			miny = min(points[static_cast<size_t>(i)].y, miny);
			maxx = max(points[static_cast<size_t>(i)].x, maxx);
			maxy = max(points[static_cast<size_t>(i)].y, maxy);
		}
		vector<double> res(4);
		res[0] = minx;
		res[1] = maxx;
		res[2] = miny;
		res[3] = maxy;
		return res;
	}

	size_t find_index(facet const& fc, int i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			if (fc.neighbors[j] == i)
				return j;
		}
		throw UniversalError("Error in find_index: Index not found");
	}
}

void Delaunay::add_point(size_t index,stack<std::pair<size_t, size_t> > &flip_stack)
{
	// Check if point is inside big triangle
	assert(InTriangle(TripleConstRef<Vector2D>(cor[olength],
		cor[olength + 1],
		cor[olength + 2]),
		cor[index]));
	const size_t triangle = Walk(index);
	const Triplet<int> outer(f[triangle].vertices);
	const Triplet<int> temp_friends(f[triangle].neighbors);
	f[triangle].vertices.set(outer.third, outer.first, static_cast<int>(index));
	f[triangle].neighbors.set(temp_friends.third,
		location_pointer + 1,
		location_pointer + 2);
	f.push_back(facet(TripleConstRef<int>(outer.first,
		outer.second,
		static_cast<int>(index)),
		TripleConstRef<int>(temp_friends.first,
			location_pointer + 2,
			static_cast<int>(triangle))));
	f.push_back(facet(TripleConstRef<int>(outer.second,
		outer.third,
		static_cast<int>(index)),
		TripleConstRef<int>(temp_friends.second,
			static_cast<int>(triangle),
			location_pointer + 1)));
	// _update the friends list of the friends
	if (temp_friends.second != last_loc)
	{
		const size_t i = find_index(f[static_cast<size_t>(temp_friends.second)], static_cast<int>(triangle));
		f[static_cast<size_t>(temp_friends.second)].neighbors[i] = location_pointer + 2;
	}
	if (temp_friends.first != last_loc)
	{
		const size_t i = find_index(f[static_cast<size_t>(temp_friends.first)], static_cast<int>(triangle));
		f[static_cast<size_t>(temp_friends.first)].neighbors[i] = location_pointer + 1;
	}
	// Calculate radius if needed
	if (CalcRadius)
	{
		radius[static_cast<size_t>(triangle)] = CalculateRadius(static_cast<int>(triangle));
		int n = int(f.size());
		int m = int(radius.size());
		if (n > m - 1)
		{
			radius.push_back(CalculateRadius(location_pointer + 1));
			radius.push_back(CalculateRadius(location_pointer + 2));
		}
		else
			if (n > m)
			{
				radius[static_cast<size_t>(location_pointer) + 1] = CalculateRadius(location_pointer + 1);
				radius.push_back(CalculateRadius(location_pointer + 2));
			}
			else
			{
				radius[static_cast<size_t>(location_pointer) + 1] = CalculateRadius(location_pointer + 1);
				radius[static_cast<size_t>(location_pointer) + 2] = CalculateRadius(location_pointer + 2);
			}
	}

	// check if flipping is needed
	flip(triangle, static_cast<size_t>(temp_friends.third),flip_stack);
	flip(static_cast<size_t>(location_pointer) + 1, static_cast<size_t>(temp_friends.first),flip_stack);
	flip(static_cast<size_t>(location_pointer) + 2, static_cast<size_t>(temp_friends.second),flip_stack);

	// _update number of facets
	location_pointer += 2;
}

void Delaunay::flip(size_t i, size_t j, stack<std::pair<size_t, size_t> > &flip_stack)
{
	if (j == static_cast<size_t>(last_loc))
		return;
	flip_stack.push(std::pair<size_t, size_t>(i, j));
	while (!flip_stack.empty())
	{
		const pair<size_t, size_t> indexes = flip_stack.top();
		// Returns the index to the point to check in coordinates and the index of the point in the facet
		const pair<int, int> check = find_diff(f[indexes.second],
			f[indexes.first]);
		const pair<int, int> other = find_diff(f[indexes.first],
			f[indexes.second]);

		facet& prefetch_1 = f[indexes.first];
		if (incircle(cor[static_cast<size_t>(prefetch_1.vertices.first)],
			cor[static_cast<size_t>(prefetch_1.vertices.second)],
			cor[static_cast<size_t>(prefetch_1.vertices.third)],
			cor[static_cast<size_t>(check.first)]) > 0)
		{
			//The point is in a circle change the facets and their friends
			const int v1 = prefetch_1.vertices[static_cast<size_t>(other.second + 1) % 3];
			const int f1 = prefetch_1.neighbors[static_cast<size_t>(other.second)];
			const int f12 = prefetch_1.neighbors[static_cast<size_t>(other.second + 2) % 3];
			facet& prefetch_2 = f[indexes.second];
			const int v2 = prefetch_2.vertices[static_cast<size_t>(check.second + 1) % 3];
			const int f2 = prefetch_2.neighbors[static_cast<size_t>(check.second + 2) % 3];
			const int f22 = prefetch_2.neighbors[static_cast<size_t>(check.second)];
			prefetch_1.vertices.set(other.first, v1, check.first);
			prefetch_2.vertices.set(check.first, v2, other.first);
			prefetch_1.neighbors.set(f1, f2, static_cast<int>(indexes.second));
			prefetch_2.neighbors.set(f22, f12, static_cast<int>(indexes.first));
			// change the friends of the friends if needed
			if (f2 != last_loc)
			{
				f[static_cast<size_t>(f2)].neighbors[static_cast<size_t>(find_index(f[static_cast<size_t>(f2)], static_cast<int>(indexes.second)))] = static_cast<int>(indexes.first);
			}
			if (f12 != last_loc)
			{
				f[static_cast<size_t>(f12)].neighbors[static_cast<size_t>(find_index(f[static_cast<size_t>(f12)], static_cast<int>(indexes.first)))] = static_cast<int>(indexes.second);
			}
			// Calculate the new radius if needed
			if (CalcRadius)
			{
				radius[indexes.first] = CalculateRadius(static_cast<int>(indexes.first));
				radius[indexes.second] = CalculateRadius(static_cast<int>(indexes.second));
			}
			// clear the checked facets
			flip_stack.pop();
			// push into the stack the new facets to check
			if (prefetch_2.neighbors.first != last_loc)
				flip_stack.push(std::pair<size_t, size_t>(indexes.second,
					static_cast<size_t>(prefetch_2.neighbors.first)));
			if (prefetch_1.neighbors.second != last_loc)
				flip_stack.push(std::pair<size_t, size_t>(indexes.first,
					static_cast<size_t>(prefetch_1.neighbors.second)));
		}
		else
		{
			// clear the checked facets
			flip_stack.pop();
		}
	}
}

void Delaunay::build_delaunay(vector<Vector2D>const& vp, vector<Vector2D> const& cpoints)
{
	cell_points = cpoints;
	DataOnlyForBuild data;
	lastFacet = 0;
	CalcRadius = false;
	length = int(vp.size() + 3);
	int len = length - 3;
	olength = static_cast<size_t>(len);
	f.clear();
	cor.clear();
	f.reserve(static_cast<size_t>(2 * length + 1 + static_cast<int>(17 * sqrt(1.*length))));
	cor.reserve(static_cast<size_t>(length + 9 * static_cast<int>(sqrt(1.*length))));
	last_loc = INT_MAX;
	for (int i = 0; i < len; i++)
	{
		cor.push_back(vp[static_cast<size_t>(i)]);
	}
	// Check point input
	CheckInput();

	// add the 3 extreme points
	Vector2D p_temp;
	vector<double> cellsize = CellSize(cell_points);
	double width = cellsize[1] - cellsize[0];
	double height = cellsize[3] - cellsize[2];
	width = max(width, height);
	height = max(width, height);
	p_temp.x = cellsize[0] - 100 * width;
	p_temp.y = cellsize[2] - 100 * height;
	cor.push_back(p_temp);
	p_temp.x = cellsize[1] + 100 * width;
	p_temp.y = cellsize[2] - 100 * height;
	cor.push_back(p_temp);
	p_temp.x = (cellsize[0] + cellsize[1]) / 2.0;
	p_temp.y = cellsize[3] + 100 * height;
	cor.push_back(p_temp);
	// Create the big triangle, and assign friends
	facet f_temp;
	f.push_back(f_temp);
	f[0].vertices[0] = len;
	f[0].vertices[1] = len + 1;
	f[0].vertices[2] = len + 2;
	for (size_t i = 0; i < 3; i++)
		f[0].neighbors[i] = last_loc;
	location_pointer = 0;
	// add the points
	size_t nloop = static_cast<size_t>(length) - 3;
	stack<std::pair<size_t, size_t> > flip_stack;
	for (size_t i = 0; i < nloop; i++)
		add_point(i,flip_stack);
	// Calculate radius
	radius.resize(f.size());
	int n = int(f.size());
	for (int i = 0; i < n; ++i)
		radius[static_cast<size_t>(i)] = CalculateRadius(i);
	CalcRadius = true;
}

double Delaunay::triangle_area(int index)
{
	const TripleConstRef<Vector2D> p
		(cor[static_cast<size_t>(f[static_cast<size_t>(index)].vertices.first)],
			cor[static_cast<size_t>(f[static_cast<size_t>(index)].vertices.second)],
			cor[static_cast<size_t>(f[static_cast<size_t>(index)].vertices.third)]);
	const double x1 = p.third.x - p.first.x;
	const double x2 = p.second.x - p.first.x;
	const double y1 = p.third.y - p.first.y;
	const double y2 = p.second.y - p.first.y;
	return -0.5*(x1*y2 - x2*y1);
}

void Delaunay::update(const vector<Vector2D>& points, vector<Vector2D>
	const& cpoints)
{
	if (logger)
		logger->output(cor, f);
	build_delaunay(points, cpoints);
}

namespace {

	int Triplet<int>::* walk_condition(const vector<Vector2D>& cor,
		const Triplet<int>& vertices,
		size_t point)
	{
		if (orient2d(TripleConstRef<Vector2D>
			(cor[static_cast<size_t>(vertices.first)],
				cor[static_cast<size_t>(vertices.second)],
				cor[point])) < 0)
			return &Triplet<int>::first;
		else if (orient2d(TripleConstRef<Vector2D>
			(cor[static_cast<size_t>(vertices.second)],
				cor[static_cast<size_t>(vertices.third)],
				cor[point])) < 0)
			return &Triplet<int>::second;
		else if (orient2d(TripleConstRef<Vector2D>
			(cor[static_cast<size_t>(vertices.third)],
				cor[static_cast<size_t>(vertices.first)],
				cor[point])) < 0)
			return &Triplet<int>::third;
		else
			return 0;
	}

	class WalkBookkeeper
	{
	public:

	private:

	};

	size_t find_new_facet(const vector<Vector2D>& cor,
		const vector<facet>& f,
		size_t point,
		size_t last_facet)
	{
		size_t res = last_facet;
		int Triplet<int>::* next = walk_condition(cor,
			f[res].vertices,
			point);
		while (next) {
			res = static_cast<size_t>(f[res].neighbors.*next);
			next = walk_condition(cor,
				f[res].vertices,
				point);
		}
		return res;
	}
}

size_t Delaunay::Walk(size_t point)
{
	lastFacet = static_cast<int>(find_new_facet(cor, f, point, static_cast<size_t>(lastFacet)));
	return static_cast<size_t>(lastFacet);
}

vector<int> Delaunay::FindContainingTetras(int StartTetra, int point)
{
	vector<int> res;
	FindContainingTetras(StartTetra, point, res);
	return res;
}

double Delaunay::FindMaxRadius(int point)
{
	const vector<int> vec = FindContainingTetras(static_cast<int>(Walk(static_cast<size_t>(point))), point);
	double r = 0;
	for (size_t i = 0; i < vec.size(); ++i)
		r = max(r, radius[static_cast<size_t>(vec[static_cast<size_t>(i)])]);
	return 2 * r;
}

void Delaunay::FindContainingTetras(int StartFacet, int point, vector<int> &result)
{
	result.clear();
	int PointLocation = FindPointInFacet(StartFacet, point);
	int NextFacet = f[static_cast<size_t>(StartFacet)].neighbors[static_cast<size_t>(PointLocation)];
	result.reserve(12);
	result.push_back(NextFacet);
	while (NextFacet != StartFacet)
	{
		PointLocation = FindPointInFacet(NextFacet, point);
		NextFacet = f[static_cast<size_t>(NextFacet)].neighbors[static_cast<size_t>(PointLocation)];
		result.push_back(NextFacet);
	}
}

int Delaunay::FindPointInFacet(int facet, int point)
{
	for (int i = 0; i < 3; ++i)
		if (f[static_cast<size_t>(facet)].vertices[static_cast<size_t>(i)] == point)
			return i;
	UniversalError eo("Error in Delaunay, FindPointInFacet");
	eo.AddEntry("Facet number", facet);
	eo.AddEntry("Point number", point);
	throw eo;
}

bool Delaunay::IsOuterFacet(int facet)const
{
	//int PointNum=length-1;
	for (int i = 0; i<3; ++i)
		for (size_t j = 0; j < 3; ++j)
			if (f[static_cast<size_t>(facet)].vertices[static_cast<size_t>(i)] == static_cast<int>(olength + j))
				return true;
	return false;
}

double Delaunay::CalculateRadius(int facet)
{
	const double big = 1e10;
	const double a = cor[static_cast<size_t>(f[static_cast<size_t>(facet)].vertices[0])].distance(cor[static_cast<size_t>(f[static_cast<size_t>(facet)].vertices[1])]);
	const double b = cor[static_cast<size_t>(f[static_cast<size_t>(facet)].vertices[0])].distance(cor[static_cast<size_t>(f[static_cast<size_t>(facet)].vertices[2])]);
	const double c = cor[static_cast<size_t>(f[static_cast<size_t>(facet)].vertices[2])].distance(cor[static_cast<size_t>(f[static_cast<size_t>(facet)].vertices[1])]);
	const double temp1 = b + c - a;
	if (temp1 <= 0)
	{
		if (a > big*b || a>big*c) // Do we have a small edge?
			return 0.5*a;
		else
			return 0.5*(b + c); // we have 3 points on a line
	}
	const double temp2 = c + a - b;
	if (temp2 <= 0)
	{
		if (b > big*a || b > big*c) // Do we have a small edge?
			return 0.5*b;
		else
			return 0.5*(a + c); // we have 3 points on a line
	}
	const double temp3 = b - c + a;
	if (temp3 <= 0)
	{
		if (c > big*b || c > big*a) // Do we have a small edge?
			return 0.5*c;
		else
			return 0.5*(b + a); // we have 3 points on a line
	}
	return a*b*c / sqrt((a + b + c)*temp1*temp2*temp3);
}

void Delaunay::CheckInput()
{
	for (size_t i = 0; i < cor.size(); ++i)
		assert(InCell(cell_points, cor[i]));
}

int Delaunay::GetOriginalIndex(int NewPoint) const
{
	return NewPoint;
}

double Delaunay::GetFacetRadius(int facet) const
{
	return radius[static_cast<size_t>(facet)];
}

void Delaunay::ChangeOlength(int n)
{
	olength = static_cast<size_t>(n);
}

void Delaunay::Changelength(int n)
{
	length = n + 3;
}

vector<Vector2D>& Delaunay::ChangeCor(void)
{
	return cor;
}

const vector<Vector2D>& Delaunay::getCor(void) const
{
	return cor;
}

const facet& Delaunay::get_facet(int index) const
{
	return f[static_cast<size_t>(index)];
}

double Delaunay::get_facet_coordinate(int Facet, int vertice, int dim)
{
	if (dim == 0)
		return cor[static_cast<size_t>(f[static_cast<size_t>(Facet)].vertices[static_cast<size_t>(vertice)])].x;
	else
		return cor[static_cast<size_t>(f[static_cast<size_t>(Facet)].vertices[static_cast<size_t>(vertice)])].y;
}

Vector2D Delaunay::get_point(size_t index) const
{
	return cor[index];
}

double Delaunay::get_cor(int index, int dim) const
{
	if (dim == 0)
		return cor[static_cast<size_t>(index)].x;
	else if (dim == 1)
		return cor[static_cast<size_t>(index)].y;
	else
		throw UniversalError("Error in Delaunay::get_cor. Invalid index");
}

int Delaunay::get_num_facet(void) const
{
	return static_cast<int>(f.size());
}

int Delaunay::get_length(void) const
{
	return length - 3;
}

int Delaunay::get_last_loc(void) const
{
	return last_loc;
}

void Delaunay::set_point(int index, Vector2D p)
{
	cor[static_cast<size_t>(index)] = p;
}

int Delaunay::GetOriginalLength(void) const
{
	return static_cast<int>(olength);
}

vector<Vector2D>& Delaunay::GetMeshPoints(void)
{
	return cor;
}

int Delaunay::GetTotalLength(void)
{
	return static_cast<int>(cor.size());
}

void Delaunay::AddBoundaryPoints(vector<Vector2D> const& points)
{
	int n = static_cast<int>(points.size());
	stack<std::pair<size_t, size_t> > flip_stack;
	//	vector<int> order=HilbertOrder(points,n);
	for (int i = 0; i < n; ++i)
	{
		cor.push_back(points[static_cast<size_t>(i)]);
		add_point(cor.size() - 1,flip_stack);
	}
}

void Delaunay::AddAditionalPoint(Vector2D const& vec)
{
	cor.push_back(vec);
}

int Delaunay::GetCorSize(void)const
{
	return static_cast<int>(cor.size());
}

bool Delaunay::IsTripleOut(int index) const
{
	int counter = 0;
	for (size_t i = 0; i < 3; ++i)
		if (IsOuterFacet(f[static_cast<size_t>(index)].neighbors[static_cast<size_t>(i)]))
			++counter;
	if (counter > 1)
		return true;
	else
		return false;
}

int Delaunay::FindTripleLoc(facet const& fct)const
{
	for (size_t i = 0; i < 3; ++i)
		if (!IsOuterFacet(fct.neighbors[static_cast<size_t>(i)]))
			return static_cast<int>((i + 1) % 3);
	throw UniversalError("Trouble in constructing boundary triangles. No inner neighbor");
}

namespace
{
	bool IsOuterQuick(facet const& f, int olength)
	{
		for (size_t i = 0; i < 3; ++i)
			if (f.vertices[static_cast<size_t>(i)] >= olength)
				return true;
		return false;
	}

	bool IsEdgeFacet(vector<facet> const& facets, facet const& f, int olength)
	{
		int counter = 0;
		for (size_t i = 0; i < 3; ++i)
		{
			if (f.vertices[static_cast<size_t>(i)] >= olength)
				return false;
			if (IsOuterQuick(facets[static_cast<size_t>(f.neighbors[static_cast<size_t>(i)])], olength))
				++counter;
		}
		if (counter > 0)
			return true;
		else
			return false;
	}

	bool CircleSegmentIntersect(Edge const& edge, Vector2D const& center, double R)
	{
		Vector2D AC = center - edge.vertices.first;
		Vector2D AB = edge.vertices.second - edge.vertices.first;
		double d = ScalarProd(AC, AB);
		if (d<0)
		{
			if (abs(AC)>R)
				return false;
			else
				return true;
		}
		double LAB = abs(AB);
		if (d > LAB*LAB)
		{
			if (abs(center - edge.vertices.second) > R)
				return false;
			else
				return true;
		}
		Vector2D closest = edge.vertices.first + AB*d / (LAB*LAB);
		if (abs(center - closest) > R)
			return false;
		else
			return true;
	}
}

vector<int> Delaunay::GetOuterFacets(int start_facet, int real_point, int olength2)
{
	int cur_facet = start_facet;
	vector<int> f_temp, containing_facets;
	f_temp.reserve(static_cast<size_t>(10 * sqrt(1.0*olength2)));
	int point_index = FindPointInFacet(cur_facet, real_point);
	if (IsOuterQuick(f[static_cast<size_t>(f[static_cast<size_t>(cur_facet)].neighbors[static_cast<size_t>(point_index)])], olength2))
	{
		point_index = (point_index + 1) % 3;
		real_point = f[static_cast<size_t>(cur_facet)].vertices[static_cast<size_t>(point_index)];
	}
	if (IsOuterQuick(f[static_cast<size_t>(f[static_cast<size_t>(cur_facet)].neighbors[static_cast<size_t>(point_index)])], olength2))
	{
		point_index = (point_index + 1) % 3;
		real_point = f[static_cast<size_t>(cur_facet)].vertices[static_cast<size_t>(point_index)];
	}
	do
	{
		FindContainingTetras(cur_facet, real_point, containing_facets);
		int old_current = cur_facet;
		for (size_t i = 0; i < containing_facets.size(); ++i)
		{
			if (IsEdgeFacet(f, f[static_cast<size_t>(containing_facets[static_cast<size_t>(i)])], olength2) &&
				containing_facets[static_cast<size_t>(i)] != old_current)
				cur_facet = containing_facets[static_cast<size_t>(i)];
			if (!IsOuterQuick(f[static_cast<size_t>(containing_facets[static_cast<size_t>(i)])], olength2))
				f_temp.push_back(containing_facets[static_cast<size_t>(i)]);
		}
		point_index = (1 + FindPointInFacet(cur_facet, real_point)) % 3;
		if (IsTripleOut(cur_facet))
			point_index = (point_index + 1) % 3;
		real_point = f[static_cast<size_t>(cur_facet)].vertices[static_cast<size_t>(point_index)];
	} while (start_facet != cur_facet);
	sort(f_temp.begin(), f_temp.end());
	f_temp = unique(f_temp);
	return f_temp;
}

vector<vector<int> > Delaunay::FindOuterPoints(vector<Edge> const& edges)
{
	// We add the points in a counter clockwise fashion
	vector<vector<int> > res(edges.size());
	if (olength < 100)
	{
		for (size_t j = 0; j < edges.size(); ++j)
		{
			res[static_cast<size_t>(j)].resize(static_cast<size_t>(olength));
			for (size_t i = 0; i < olength; ++i)
				res[static_cast<size_t>(j)][i] = static_cast<int>(i);
		}
		return res;
	}
	vector<int> res_temp, outer_points, f_temp, f_add(f.size(), 0);
	res_temp.reserve(static_cast<size_t>(20 * sqrt(1.0*static_cast<double>(olength))));
	f_temp.reserve(static_cast<size_t>(10 * sqrt(1.0*static_cast<double>(olength))));
	outer_points.reserve(static_cast<size_t>(10 * sqrt(1.0*static_cast<double>(olength))));
	// Walk to an outer point
	int cur_facet = static_cast<int>(Walk(olength));
	vector<vector<int> > toduplicate(edges.size());
	vector<bool> checked(f.size(), false);
	AddOuterFacets(cur_facet, toduplicate, edges, checked);
	for (size_t i = 0; i < edges.size(); ++i)
	{
		sort(toduplicate[static_cast<size_t>(i)].begin(), toduplicate[static_cast<size_t>(i)].end());
		toduplicate[static_cast<size_t>(i)] = unique(toduplicate[static_cast<size_t>(i)]);
	}
	return toduplicate;
}

void Delaunay::AddRigid(vector<Edge> const& edges,
	vector<vector<int> > &toduplicate)
{
	vector<int> toremove;
	for (size_t i = 0; i < edges.size(); ++i)
	{
		toremove.clear();
		if (toduplicate[i].empty())
			continue;
		vector<Vector2D> toadd;
		toadd.reserve(toduplicate[i].size());
		Vector2D par(Parallel(edges[i]));
		par = par / abs(par);
		for (size_t j = 0; j < toduplicate[i].size(); ++j)
		{
			Vector2D temp = cor[static_cast<size_t>(toduplicate[i][j])] - edges[i].vertices.first;
			temp = 2 * par*ScalarProd(par, temp) - temp + edges[i].vertices.first;
			if (InTriangle(TripleConstRef<Vector2D>
				(cor[static_cast<size_t>(olength)],
					cor[static_cast<size_t>(olength + 1)],
					cor[static_cast<size_t>(olength + 2)]),
				temp))
				toadd.push_back(temp);
			else
				toremove.push_back(static_cast<int>(j));
		}
		RemoveVector(toduplicate[i], toremove);
		vector<int> order = HilbertOrder(toadd, static_cast<int>(toadd.size()));
		ReArrangeVector(toadd, order);
		try
		{
			AddBoundaryPoints(toadd);
		}
		catch (UniversalError &eo)
		{
			eo.AddEntry("Error in AddRigid", 0);
			throw;
		}
		ReArrangeVector(toduplicate[i], order);
	}
}

namespace
{
	vector<Edge> GetCornerEdges(OuterBoundary const* obc)
	{
		const double dx = obc->GetGridBoundary(Right) - obc->GetGridBoundary(Left);
		const double dy = obc->GetGridBoundary(Up) - obc->GetGridBoundary(Down);
		vector<Edge> res;
		const Vector2D RU(obc->GetGridBoundary(Right), obc->GetGridBoundary(Up));
		const Vector2D LU(obc->GetGridBoundary(Left), obc->GetGridBoundary(Up));
		const Vector2D LD(obc->GetGridBoundary(Left), obc->GetGridBoundary(Down));
		const Vector2D RD(obc->GetGridBoundary(Right), obc->GetGridBoundary(Down));
		res.push_back(Edge(RU, Vector2D(dx, 0) + RU, 0, 0));
		res.push_back(Edge(RU, Vector2D(0, dy) + RU, 0, 0));
		res.push_back(Edge(LU, Vector2D(0, dy) + LU, 0, 0));
		res.push_back(Edge(LU, Vector2D(-dx, 0) + LU, 0, 0));
		res.push_back(Edge(LD, Vector2D(-dx, 0) + LD, 0, 0));
		res.push_back(Edge(LD, Vector2D(0, -dy) + LD, 0, 0));
		res.push_back(Edge(RD, Vector2D(0, -dy) + RD, 0, 0));
		res.push_back(Edge(RD, Vector2D(dx, 0) + RD, 0, 0));
		return res;
	}
}

vector<vector<int> > Delaunay::AddPeriodic(OuterBoundary const* obc, vector<Edge> const& edges,
	vector<vector<int> > &toduplicate)
{
	const double dx = obc->GetGridBoundary(Right) - obc->GetGridBoundary(Left);
	const double dy = obc->GetGridBoundary(Up) - obc->GetGridBoundary(Down);
	for (size_t i = 0; i < edges.size(); ++i)
	{
		if (toduplicate[static_cast<size_t>(i)].empty())
			continue;
		Vector2D change;
		switch (i)
		{
		case(0) :
			change.x = -dx;
			break;
		case(1) :
			change.y = -dy;
			break;
		case(2) :
			change.x = dx;
			break;
		case(3) :
			change.y = dy;
			break;
		}
		vector<Vector2D> toadd;
		toadd.reserve(toduplicate[static_cast<size_t>(i)].size());
		//vector<int> pointstemp(toduplicate[static_cast<size_t>(i)].size());
		for (size_t j = 0; j < toduplicate[static_cast<size_t>(i)].size(); ++j)
		{
			toadd.push_back(cor[static_cast<size_t>(toduplicate[static_cast<size_t>(i)][static_cast<size_t>(j)])] + change);
			//	pointstemp[j]=j;
		}
		vector<int> order = HilbertOrder(toadd, static_cast<int>(toadd.size()));
		ReArrangeVector(toadd, order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(toduplicate[static_cast<size_t>(i)], order);
		//toduplicate[i]=pointstemp;
	}
	// Done with sides do corners now
	vector<Edge> corneredges = GetCornerEdges(obc);
	vector<vector<int> > corners(toduplicate.size());
	for (size_t i = 0; i < toduplicate.size(); ++i)
	{
		for (size_t j = 0; j < toduplicate[static_cast<size_t>(i)].size(); ++j)
		{
			const int facet_loc = static_cast<int>(Walk(static_cast<size_t>(toduplicate[static_cast<size_t>(i)][static_cast<size_t>(j)])));
			const Vector2D center = cor[static_cast<size_t>(toduplicate[static_cast<size_t>(i)][static_cast<size_t>(j)])];
			const double R = 2 * GetMaxRadius(toduplicate[static_cast<size_t>(i)][static_cast<size_t>(j)], facet_loc);
			if (CircleSegmentIntersect(corneredges[2 * static_cast<size_t>(i)], center, R))
				corners[static_cast<size_t>(i)].push_back(toduplicate[static_cast<size_t>(i)][static_cast<size_t>(j)]);
			if (CircleSegmentIntersect(corneredges[(2 * static_cast<size_t>(i) + 7) % 8], center, R))
				corners[(static_cast<size_t>(i) + 3) % 4].push_back(toduplicate[static_cast<size_t>(i)][static_cast<size_t>(j)]);
		}
	}
	for (size_t i = 0; i < corners.size(); ++i)
	{
		if (corners[static_cast<size_t>(i)].empty())
			continue;
		sort(corners[static_cast<size_t>(i)].begin(), corners[static_cast<size_t>(i)].end());
		corners[static_cast<size_t>(i)] = unique(corners[static_cast<size_t>(i)]);
		Vector2D change;
		switch (i)
		{
		case(0) :
			change.x = -dx;
			change.y = -dy;
			break;
		case(1) :
			change.y = -dy;
			change.x = dx;
			break;
		case(2) :
			change.x = dx;
			change.y = dy;
			break;
		case(3) :
			change.y = dy;
			change.x = -dx;
			break;
		}
		vector<Vector2D> toadd;
		toadd.reserve(corners[static_cast<size_t>(i)].size());
		//		vector<int> pointstemp(corners[i].size());
		for (size_t j = 0; j < corners[static_cast<size_t>(i)].size(); ++j)
		{
			toadd.push_back(cor[static_cast<size_t>(corners[static_cast<size_t>(i)][static_cast<size_t>(j)])] + change);
			//		pointstemp[j]=j;
		}
		vector<int> order = HilbertOrder(toadd, static_cast<int>(toadd.size()));
		ReArrangeVector(toadd, order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(corners[static_cast<size_t>(i)], order);
		//	corners[i]=pointstemp;
	}
	return corners;
}

void Delaunay::AddHalfPeriodic(OuterBoundary const* obc, vector<Edge> const& edges,
	vector<vector<int> > &toduplicate)
{
	const double dx = obc->GetGridBoundary(Right) - obc->GetGridBoundary(Left);
	//	const double dy=obc->GetGridBoundary(Up)-obc->GetGridBoundary(Down);
	for (size_t i = 0; i < edges.size(); ++i)
	{
		if (toduplicate[static_cast<size_t>(i)].empty())
			continue;
		Vector2D change;
		switch (i)
		{
		case(0) :
			change.x = -dx;
			break;
		case(1) :
			break;
		case(2) :
			change.x = dx;
			break;
		case(3) :
			break;
		}
		vector<Vector2D> toadd;
		toadd.reserve(toduplicate[static_cast<size_t>(i)].size());
		//vector<int> pointstemp(toduplicate[static_cast<size_t>(i)].size());
		Vector2D par(Parallel(edges[static_cast<size_t>(i)]));
		par = par / abs(par);
		for (size_t j = 0; j < toduplicate[static_cast<size_t>(i)].size(); ++j)
		{
			Vector2D temp = cor[static_cast<size_t>(toduplicate[static_cast<size_t>(i)][static_cast<size_t>(j)])];
			if (i % 2 == 1)
			{
				temp -= edges[static_cast<size_t>(i)].vertices.first;
				temp = 2 * par*ScalarProd(par, temp) - temp + edges[static_cast<size_t>(i)].vertices.first;
			}
			toadd.push_back(temp + change);
			//pointstemp[j]=j;
		}
		vector<int> order = HilbertOrder(toadd, static_cast<int>(toadd.size()));
		ReArrangeVector(toadd, order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(toduplicate[static_cast<size_t>(i)], order);
		//toduplicate[i]=pointstemp;
	}
}

vector<vector<int> > Delaunay::BuildBoundary(OuterBoundary const* obc, vector<Edge> const& edges)
{
	vector<vector<int> > toduplicate = FindOuterPoints(edges);
#ifdef RICH_MPI
	OrgIndex.clear();
#endif
	if (obc->GetBoundaryType() == Rectengular)
	{
		AddRigid(edges, toduplicate);
#ifdef RICH_MPI
		for (size_t i = 0; i < toduplicate.size(); ++i)
			for (size_t j = 0; j < toduplicate[i].size(); ++j)
				OrgIndex.push_back(toduplicate[i][j]);
#endif
	}
	else
	{
		if (obc->GetBoundaryType() == Periodic)
		{
			vector<vector<int> > corners = AddPeriodic(obc, edges, toduplicate);
			for (size_t i = 0; i < 4; ++i)
				toduplicate.push_back(corners[static_cast<size_t>(i)]);
		}
		else
		{
			AddHalfPeriodic(obc, edges, toduplicate);
		}
	}
	return toduplicate;
}

Vector2D Delaunay::GetCircleCenter(int index)const
{
	Vector2D center;
	facet const& F = f[static_cast<size_t>(index)];
	double x1 = cor[static_cast<size_t>(F.vertices[0])].x;
	double x2 = cor[static_cast<size_t>(F.vertices[1])].x;
	double x3 = cor[static_cast<size_t>(F.vertices[2])].x;
	double y1 = cor[static_cast<size_t>(F.vertices[0])].y;
	double y2 = cor[static_cast<size_t>(F.vertices[1])].y;
	double y3 = cor[static_cast<size_t>(F.vertices[2])].y;
	// Do we have a case where two point are very close compared to the third?
	double d12 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
	double d23 = (x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2);
	double d13 = (x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3);
	int scenario = 0;
	if (d12 < 0.1*(d23 + d13))
		scenario = 1;
	else
		if (d23 < 0.1*(d13 + d12))
			scenario = 3;
		else
			if (d13 < 0.1*(d23 + d12))
				scenario = 2;
	switch (scenario)
	{
	case(0) :
	case(1) :
	case(2) :
	{
		x2 -= x1;
		x3 -= x1;
		y2 -= y1;
		y3 -= y1;
		double d_inv = 1 / (2 * (x2*y3 - y2*x3));
		center.Set((y3*(x2*x2 + y2*y2) - y2*(x3*x3 + y3*y3))*d_inv + x1,
			(-x3*(x2*x2 + y2*y2) + x2*(x3*x3 + y3*y3))*d_inv + y1);
		break;
	}
	case(3) :
	{
		x1 -= x2;
		x3 -= x2;
		y1 -= y2;
		y3 -= y2;
		double d_inv = 1 / (2 * (x3*y1 - y3*x1));
		center.Set((y1*(x3*x3 + y3*y3) - y3*(x1*x1 + y1*y1))*d_inv + x2,
			(x3*(x1*x1 + y1*y1) - x1*(x3*x3 + y3*y3))*d_inv + y2);
		break;
	}
	default:
		throw UniversalError("Unhandled case in switch statement VoronoiMesh::get_center");
	}
	return center;
}

double Delaunay::GetMaxRadius(int point, int startfacet)
{
	double res = 0;
	vector<int> neigh = FindContainingTetras(startfacet, point);
	for (size_t i = 0; i < neigh.size(); ++i)
		res = max(res, radius[static_cast<size_t>(neigh[static_cast<size_t>(i)])]);
	return res;
}

void Delaunay::AddOuterFacets(int tri, vector<vector<int> > &toduplicate,
	vector<Edge> const& edges, vector<bool> &checked)
{
	stack<int> tocheck;
	tocheck.push(tri);
	while (!tocheck.empty())
	{
		int cur_facet = tocheck.top();
		tocheck.pop();
		for (size_t i = 0; i < 3; ++i)
		{
			bool added = false;
			if (checked[static_cast<size_t>(f[static_cast<size_t>(cur_facet)].vertices[static_cast<size_t>(i)])] || (f[static_cast<size_t>(cur_facet)].vertices[static_cast<size_t>(i)] >= static_cast<int>(olength)))
				continue;
			vector<int> neigh = FindContainingTetras(cur_facet, f[static_cast<size_t>(cur_facet)].vertices[static_cast<size_t>(i)]);
			for (size_t k = 0; k < neigh.size(); ++k)
			{
				Vector2D center = GetCircleCenter(neigh[static_cast<size_t>(k)]);
				for (size_t l = 0; l < edges.size(); ++l)
				{
					if (CircleSegmentIntersect(edges[static_cast<size_t>(l)], center, radius[static_cast<size_t>(neigh[static_cast<size_t>(k)])]))
					{
						toduplicate[static_cast<size_t>(l)].push_back(f[static_cast<size_t>(cur_facet)].vertices[static_cast<size_t>(i)]);
						added = true;
					}
				}
			}
			checked[static_cast<size_t>(f[static_cast<size_t>(cur_facet)].vertices[static_cast<size_t>(i)])] = true;
			if (added)
			{
				for (size_t j = 0; j < neigh.size(); ++j)
					tocheck.push(neigh[static_cast<size_t>(j)]);
			}
		}
	}
}

#ifdef RICH_MPI
int Delaunay::findSomeOuterPoint(void)
{
	const size_t cur_facet = Walk(olength);
	for (size_t i = 0; i < 3; ++i) {
		const int candidate = f.at(cur_facet).vertices[i];
		if (candidate < static_cast<int>(olength))
			return candidate;
	}
	assert(false && "something went wrong");
}

namespace 
{
	vector<int>
		calc_neighbors_own_edges
		(const Tessellation& t_proc,
			const vector<Edge>& edge_list)
	{
		vector<int> res;
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		BOOST_FOREACH(const Edge& edge, edge_list) 
		{
			const int other = (edge.neighbors.first + edge.neighbors.second) -	rank;
			if (other < t_proc.GetPointNo()) {
				res.push_back(other);
			}
		}
		return res;
	}

	stack<int> initialise_tocheck
		(const vector<int>& neightemp)
	{
		stack<int> res;
		for (size_t i = 0; i < neightemp.size(); ++i)
			res.push(neightemp[i]);
		return res;
	}

	vector<int> calc_self_intersection
		(const vector<Edge>& edge_list,
			const Circle& circle)
	{
		vector<int> res;
		for (size_t i = 0; i < edge_list.size(); ++i) {
			const Edge& edge = edge_list.at(i);
			if (edge_circle_intersect(edge, circle))
				res.push_back(static_cast<int>(i));
		}
		return res;
	}
}

vector<vector<int> > Delaunay::AddOuterFacetsMPI
(int point,
	vector<vector<int> > &toduplicate,
	vector<int> &neigh,
	vector<bool> &checked,
	Tessellation const &tproc,
	const vector<Edge>& own_edges,
	bool recursive)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<vector<int> > res;
	vector<int> vtemp;
	if (!recursive)
		res.resize(own_edges.size());
	stack<int> tocheck = initialise_tocheck
		(FindContainingTetras
			(static_cast<int>(Walk(static_cast<size_t>(point))), point));
	if (recursive)
	{
		vector<int> allouter;
		for (size_t i = 0; i < toduplicate.size(); ++i)
		{
			for (size_t j = 0; j < toduplicate[i].size(); ++j)
			{
				vector<int> temp = FindContainingTetras
					(static_cast<int>(Walk(static_cast<size_t>(toduplicate[i][j]))), toduplicate[i][j]);
				for (size_t k = 0; k < temp.size(); ++k)
					allouter.push_back(temp[k]);
			}
		}
		sort(allouter.begin(), allouter.end());
		allouter = unique(allouter);
		for (size_t i = 0; i < allouter.size(); ++i)
			tocheck.push(allouter[i]);
	}
	while (!tocheck.empty())
	{
		int cur_facet = tocheck.top();
		tocheck.pop();
		for (size_t i = 0; i < 3; ++i)
		{
			bool added = false;
			int max_neigh = 0;
			if (f[static_cast<size_t>(cur_facet)].vertices[i] >=
				static_cast<int>(olength))
				continue;
			if (checked[static_cast<size_t>
				(f[static_cast<size_t>(cur_facet)].vertices[i])])
				continue;
			vector<int> neighs = FindContainingTetras(cur_facet, f[static_cast<size_t>(cur_facet)].vertices[i]);
			for (size_t k = 0; k < neighs.size(); ++k)
			{
				Circle circ(GetCircleCenter(neighs[k]), radius[static_cast<size_t>(neighs[k])]);
				vector<int> cputosendto;
				if (recursive)
					find_affected_cells_recursive(tproc,rank, circ, cputosendto);
				else
					cputosendto = find_affected_cells
					(tproc, rank, circ,vtemp);
				sort(cputosendto.begin(), cputosendto.end());
				cputosendto = unique(cputosendto);

				RemoveVal(cputosendto,rank);
				if (!recursive) 
				{
					const vector<int> self_intersection = calc_self_intersection(own_edges, circ);
					if (!self_intersection.empty())
						added = true;
					BOOST_FOREACH(int sindex, self_intersection)
						res[static_cast<size_t>(sindex)].push_back(f[static_cast<size_t>(cur_facet)].vertices[i]);
				}
				else
				{
					for (size_t jj = 0; jj < 3; ++jj)
						max_neigh = max(max_neigh, f[static_cast<size_t>(neighs[k])].vertices[jj]);
				}
				if (!cputosendto.empty())
				{
					added = true;
					for (size_t j = 0; j < cputosendto.size(); ++j)
					{
						size_t index = static_cast<size_t>(find(neigh.begin(), neigh.end(), cputosendto[j])
							- neigh.begin());
						if (index < neigh.size())
							toduplicate.at(index).push_back(f[static_cast<size_t>(cur_facet)].vertices[i]);
						else {
							neigh.push_back(cputosendto[j]);
							toduplicate.push_back
								(vector<int>(1, f[static_cast<size_t>(cur_facet)].vertices[i]));
						}
					}
				}
			}
			checked[static_cast<size_t>(f[static_cast<size_t>(cur_facet)].vertices[i])] = true;
			if (added||(recursive&&static_cast<size_t>(max_neigh)>=olength))
			{
				for (size_t j = 0; j < neighs.size(); ++j)
				{
					//if (!IsOuterFacet(neighs[j]))
						tocheck.push(neighs[j]);
				}
			}
		}
	}
	return res;
}

pair<vector<vector<int> >, vector<vector<int> > >
Delaunay::findOuterPoints
(const Tessellation& t_proc,
	const vector<Edge>& edge_list,
	const vector<Edge>& box_edges,
	vector<vector<int> > &NghostIndex)
{
	vector<int> neighbors_own_edges =
		calc_neighbors_own_edges(t_proc, edge_list);
	const size_t some_outer_point = findSomeOuterPoint();

	vector<vector<int> > to_duplicate(neighbors_own_edges.size());
	vector<bool> checked(olength, false);
	vector<vector<int> > self_points =
		AddOuterFacetsMPI
		(static_cast<int>(some_outer_point),
			to_duplicate, // indices of points to send
			neighbors_own_edges, // Rank of processes to send to
			checked,
			t_proc,
			box_edges);
	BOOST_FOREACH(vector<int>& line, to_duplicate) 
	{
		sort(line.begin(), line.end());
		line = unique(line);
	}

	// Communication
	vector<vector<Vector2D> > incoming(neighbors_own_edges.size());
	vector<MPI_Request> req(neighbors_own_edges.size());
	vector<vector<double> > tosend(neighbors_own_edges.size());
	double dtemp = 0;
	for (size_t i = 0; i < neighbors_own_edges.size(); ++i)
	{
		const int dest = neighbors_own_edges.at(i);
		tosend[i] = list_serialize(VectorValues(cor, to_duplicate[i]));
		int size = static_cast<int>(tosend[i].size());
		if (size == 0)
			MPI_Isend(&dtemp, 1, MPI_DOUBLE,dest, 1, MPI_COMM_WORLD, &req[i]);
		else
		{
			if (size < 2)
				throw UniversalError("Wrong send size");
			MPI_Isend(&tosend[i][0], size, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &req[i]);
		}
	}
	for (size_t i = 0; i < neighbors_own_edges.size(); ++i)
	{
		vector<double> temprecv;
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		temprecv.resize(static_cast<size_t>(max(count, 1)));
		MPI_Recv(&temprecv[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (status.MPI_TAG == 0)
		{
			size_t location = static_cast<size_t>(std::find(neighbors_own_edges.begin(), neighbors_own_edges.end(),
				status.MPI_SOURCE) - neighbors_own_edges.begin());
			if (location >= neighbors_own_edges.size())
				throw UniversalError("Bad location in mpi exchange");
			try
			{
				incoming[location] = list_unserialize(temprecv, cor[0]);
			}
			catch (UniversalError &eo)
			{
				eo.AddEntry("Error in first send in triangulation", 0.0);
				throw;
			}
		}
		else
			if (status.MPI_TAG != 1)
				throw UniversalError("Wrong mpi tag");
	}
	MPI_Waitall(static_cast<int>(req.size()), &req[0], MPI_STATUSES_IGNORE);
	// Incorporate points recieved into triangulation
	BOOST_FOREACH(vector<int> &line, self_points)
	{
		sort(line.begin(), line.end());
		line = unique(line);
	}
	AddRigid(box_edges,self_points);
	for (size_t i = 0; i < self_points.size(); ++i)
		for (size_t j = 0; j < self_points[i].size(); ++j)
			OrgIndex.push_back(self_points[i][j]);

	NghostIndex.clear();
	NghostIndex.resize(incoming.size());
	for (size_t i = 0; i < incoming.size(); ++i)
	{
		for (size_t j = 0; j < incoming.at(i).size(); ++j)
		{
			OrgIndex.push_back(static_cast<int>(cor.size() + j));
			NghostIndex[i].push_back(static_cast<int>(cor.size() + j));
		}
		AddBoundaryPoints(incoming.at(i));
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return pair<vector<vector<int> >, vector<vector<int> > >
		(to_duplicate, self_points);
}

namespace {
	template<class T> bool is_in
		(const T& t,
			const vector<T>& v)
	{
		BOOST_FOREACH(const T&m, v) {
			if (t == m)
				return true;
		}
		return false;
	}
}

vector<vector<int> >
Delaunay::boundary_intersection_check
(const vector<Edge>& edges,
	const vector<vector<int> >& to_duplicate)
{
	vector<vector<int> > res;
	res.reserve(edges.size());
	BOOST_FOREACH(const Edge& edge, edges) {
		res.push_back(vector<int>());
		BOOST_FOREACH(const vector<int>& line, to_duplicate) {
			BOOST_FOREACH(const int index, line) {
				const Circle circle
					(GetCircleCenter(index), radius[index]);
				if (edge_circle_intersect(edge, circle))
					res.back().push_back(index);
			}
		}
	}
	return res;
}

pair<vector<vector<int> >, vector<int> > Delaunay::FindOuterPoints2
(const Tessellation& t_proc,
	const vector<Edge>& edge_list,
	vector<vector<int> > &to_duplicate,
	vector<vector<int> >& self_points,
	const vector<Edge>& box_edges,
	vector<vector<int> > &NghostIndex)
{
	const vector<vector<int> > boundary_points =
		boundary_intersection_check(box_edges, to_duplicate);
	vector<vector<int> > real_boundary_points
		(box_edges.size());
	for (size_t i = 0; i < boundary_points.size(); ++i)
	{
		sort(self_points.at(i).begin(), self_points.at(i).end());
		BOOST_FOREACH(int bp, boundary_points.at(i))
		{
			if (!binary_search(self_points.at(i).begin(),
				self_points.at(i).end(),
				bp))
			{
				real_boundary_points.at(i).push_back(bp);
			}
		}
		sort(real_boundary_points.at(i).begin(),
			real_boundary_points.at(i).end());
		real_boundary_points.at(i) = unique(real_boundary_points.at(i));
	}

	vector<vector<int> > to_duplicate_2 = to_duplicate;
	BOOST_FOREACH(vector<int>& line, to_duplicate_2)
		sort(line.begin(), line.end());
	vector<int> neighbors_own_edges =
		calc_neighbors_own_edges(t_proc, edge_list);
	vector<int> old_neighbors = neighbors_own_edges;
	assert(!to_duplicate.empty());
	BOOST_FOREACH(const vector<int>& line, to_duplicate)
		assert(!line.empty());
	vector<bool> checked(olength, false);
	const size_t some_outer_point = to_duplicate[0][0];
	AddOuterFacetsMPI
		(static_cast<int>(some_outer_point),
			to_duplicate, // indices of points to send
			neighbors_own_edges, // Rank of processes to send to
			checked,
			t_proc,
			edge_list,
			true); // recursive

	// Communication
	int wsize;
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	vector<int> totalk(static_cast<size_t>(wsize), 0);
	vector<int> scounts(totalk.size(), 1);
	for (size_t i = 0; i < neighbors_own_edges.size(); ++i)
		totalk[neighbors_own_edges[i]] = 1;
	int nrecv;
	MPI_Reduce_scatter(&totalk[0], &nrecv, &scounts[0], MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);

	vector<MPI_Request> req(neighbors_own_edges.size());
	for (size_t i = 0; i < neighbors_own_edges.size(); ++i)
		MPI_Isend(&wsize, 1, MPI_INT, neighbors_own_edges[i], 3, MPI_COMM_WORLD, &req[i]);
	vector<int> talkwithme;
	for (int i = 0; i < nrecv; ++i)
	{
		MPI_Status status;
		MPI_Recv(&wsize, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
		talkwithme.push_back(status.MPI_SOURCE);
	}
	vector<size_t> indices;
	for (size_t i = 0; i < neighbors_own_edges.size(); ++i)
	{
		if (is_in(neighbors_own_edges[i],talkwithme))
			indices.push_back(i);
	}

	// Symmetrisation
	neighbors_own_edges =
		VectorValues
		(neighbors_own_edges,
			indices);
	to_duplicate =
		VectorValues
		(to_duplicate,
			indices);
	// Get rid of duplicate points
	for (size_t i = 0; i < to_duplicate.size(); ++i)
	{
		sort(to_duplicate[i].begin(), to_duplicate[i].end());
		to_duplicate[i] = unique(to_duplicate[i]);
	}
	vector<vector<int> > messages(to_duplicate.size());
	for (size_t i = 0; i < to_duplicate.size(); ++i) {
		const vector<int>::const_iterator it =
			find(old_neighbors.begin(),
				old_neighbors.end(),
				neighbors_own_edges.at(i));
		for (size_t j = 0; j < to_duplicate.at(i).size(); ++j) {

			if (it != old_neighbors.end()) {
				const size_t my_index = static_cast<size_t>(it - old_neighbors.begin());
				if (!binary_search
					(to_duplicate_2.at(my_index).begin(),
						to_duplicate_2.at(my_index).end(),
						to_duplicate.at(i).at(j)))
					messages.at(i).push_back(to_duplicate.at(i).at(j));
			}
			else
				messages.at(i).push_back(to_duplicate.at(i).at(j));
		}
	}
	MPI_Waitall(static_cast<int>(neighbors_own_edges.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	// Point exchange
	req.clear();
	req.resize(neighbors_own_edges.size());
	double dtemp = 0;
	vector<vector<Vector2D> > incoming(neighbors_own_edges.size());
	vector<vector<double> > tosend(neighbors_own_edges.size());
	for (size_t i = 0; i < neighbors_own_edges.size(); ++i)
	{
		const int dest = neighbors_own_edges.at(i);
		tosend[i] = list_serialize(VectorValues(cor, messages.at(i)));
		int size = static_cast<int>(tosend[i].size());
		if (size > 0)
		{
			if (size < 2)
				throw UniversalError("Wrong send size");
			MPI_Isend(&tosend[i][0], size, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD, &req[i]);
		}
		else
			MPI_Isend(&dtemp, 1, MPI_DOUBLE, dest,1, MPI_COMM_WORLD, &req[i]);
	}
	for (size_t i = 0; i < neighbors_own_edges.size(); ++i)
	{
		vector<double> temprecv;
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		temprecv.resize(static_cast<size_t>(count));
		MPI_Recv(&temprecv[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (status.MPI_TAG == 0)
		{
			size_t location = static_cast<size_t>(std::find(neighbors_own_edges.begin(), neighbors_own_edges.end(), status.MPI_SOURCE) -
				neighbors_own_edges.begin());
			if (location >= neighbors_own_edges.size())
				throw UniversalError("Bad location in mpi exchange");
			try
			{
				incoming[location] = list_unserialize(temprecv, cor[0]);
			}
			catch (UniversalError &eo)
			{
				eo.AddEntry("Error in second send in triangulation", 0.0);
				eo.AddEntry("Mpi status", static_cast<double>(status.MPI_SOURCE));
				eo.AddEntry("Mpi tag", static_cast<double>(status.MPI_TAG));
				eo.AddEntry("Mpi count", static_cast<double>(count));
				throw;
			}
		}
		else
			if (status.MPI_TAG != 1)
				throw UniversalError("Wrong mpi tag");
	}
	MPI_Waitall(static_cast<int>(req.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	// Incorporate points recieved into triangulation
	NghostIndex.resize(incoming.size());
	for (size_t i = 0; i < incoming.size(); ++i)
	{
		for (size_t j = 0; j < incoming.at(i).size(); ++j)
		{
			NghostIndex[i].push_back(static_cast<int>(cor.size() + j));
			OrgIndex.push_back(static_cast<int>(cor.size() + j));
		}
		AddBoundaryPoints(incoming.at(i));
	}

	AddRigid(box_edges,real_boundary_points);
	for (size_t i = 0; i < real_boundary_points.size(); ++i)
		for (size_t j = 0; j < real_boundary_points[i].size(); ++j)
			OrgIndex.push_back(real_boundary_points[i][j]);

	for (size_t i = 0; i < self_points.size(); ++i)
		if (!real_boundary_points[i].empty())
			self_points[i].insert(self_points[i].end(), real_boundary_points[i].begin(),
				real_boundary_points[i].end());

	return  pair<vector<vector<int> >,vector<int> > (to_duplicate,neighbors_own_edges);
}

pair<vector<vector<int> >, vector<int> > Delaunay::BuildBoundary
(OuterBoundary const* obc,
	Tessellation const& tproc,
	vector<vector<int> >& Nghost)
{
	vector<Edge> edges;
	OrgIndex.clear();
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<int> edge_index = tproc.GetCellEdges(rank);
	for (size_t i = 0; i < edge_index.size(); ++i)
		edges.push_back(tproc.GetEdge(edge_index[i]));
	vector<Edge> box_edges = obc->GetBoxEdges();
	pair<vector<vector<int> >, vector<vector<int> > > to_duplicate =
		findOuterPoints(tproc, edges, box_edges, Nghost);
	return FindOuterPoints2(tproc,edges,to_duplicate.first, to_duplicate.second,box_edges, Nghost);
}

int Delaunay::GetOrgIndex(int index)const
{
	if (index < static_cast<int>(olength))
		return static_cast<int>(olength);
	else
		return OrgIndex.at(index - 3 - static_cast<int>(olength));
}

#endif // RICH_MPI
