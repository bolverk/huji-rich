#include "Delaunay.hpp"
#include <vector>
#include <cmath>
#include "../misc/triplet.hpp"
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <numeric>
#include <boost/range/combine.hpp>
#ifdef RICH_MPI
#include <mpi.h>
#endif // RICH_MPI

namespace {

  template<class T> bool is_in(const T& t, const Triplet<T>& trip)
  {
    for(auto s : trip){
      if(t==s)
	return true;
    }
    return false;
  }

  pair<int, int> find_diff(const facet& f1, const facet& f2)
  {
    for(size_t i=0;i<3;++i){
      if(not is_in(f1.vertices[i],f2.vertices))
	return pair<int,int>(f1.vertices[i],i);
    }
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
	logger(nullptr)
#ifdef RICH_MPI
	,OrgIndex(vector<size_t>())
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
	bool InCell(const vector<Vector2D>& points, 
		    const Vector2D& p)
	{
	  for (size_t i = 0, n = points.size(); i < n; ++i){
	    if (CrossProduct(points[i] - p, points[(i + 1) % n] - p) < 0)
	      return false;
	  }
	  return true;
	}

  template<class S, class T> T convert_reduce
  (const vector<S>& slist,
   const std::function<T(S)>& f_con,
   const std::function<T(T,T)>& f_red)
  {
    T res = f_con(slist.front());
    for(size_t i=1;i<slist.size();++i)
      res = f_red(res, f_con(slist[i]));
    return res;
  }

  std::array<double,4> CellSize(vector<Vector2D> const& points)
	{
	  std::array<double,4> res;
	  size_t const n = points.size();
	  double minx = points[0].x;
	  double miny = points[0].y;
	  double maxx = minx;
	  double maxy = miny;
	  for (size_t i = 1; i < n; ++i)
	  {
		  minx = std::min(points[i].x, minx);
		  miny = std::min(points[i].y, miny);
		  maxx = std::max(points[i].x, maxx);
		  maxy = std::max(points[i].y, maxy);
	  }
	  res[0] = minx;
	  res[1] = maxx;
	  res[2] = miny;
	  res[3] = maxy;
	  return res;
	}

	size_t find_index(facet const& fc, size_t i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
		  if (fc.neighbors[j] == i)
				return j;
		}
		throw UniversalError("Error in find_index: Index not found");
	}

  template<class T1, class T2, int N> std::array<pair<T1, T2>, N> zip2
  (const std::array<T1, N>& v1, 
   const std::array<T2, N>& v2)
  {
    std::array<pair<T1, T2>, N> res;
    transform(v1.begin(),
	      v1.end(),
	      v2.begin(),
	      res.begin(),
	      [](const T1& t1, const T2& t2)
	      {return pair<T1, T2>(t1, t2);});
    return res;
  }
}

void Delaunay::check_if_flipping_is_needed
(size_t triangle,
 const Triplet<size_t>& temp_friends,
 stack<std::pair<size_t, size_t> >& flip_stack)
{
  for(const auto& zipped : 
	zip2<size_t,size_t,3>
	({triangle, 
	    location_pointer+1,
	    location_pointer+2},
	{temp_friends.third, 
	    temp_friends.first,
	    temp_friends.second}))
    flip(zipped.first, zipped.second, flip_stack);
}

void Delaunay::update_friends_of_friends
(size_t triangle, const Triplet<size_t>& temp_friends)
{
  for(const auto& itm :
    {pair<size_t,size_t>(temp_friends.first, 1),
	pair<size_t,size_t>(temp_friends.second, 2)})
    if(itm.first != last_loc){
      const size_t i = find_index
	(f[itm.first],
	 triangle);
      f[itm.first].neighbors[i] = location_pointer + itm.second;
    }
}

void Delaunay::update_f_in_add_point
(size_t triangle, const Triplet<size_t>& temp_friends, size_t index)
{
  const Triplet<size_t> outer(f[triangle].vertices);
  f[triangle].vertices.set
    (outer.third,
     outer.first,
     index);
  f[triangle].neighbors.set(temp_friends.third,
			    location_pointer + 1,
			    location_pointer + 2);
  const facet facet1
    (TripleConstRef<size_t>
     (outer.first,
      outer.second,
      index),
     TripleConstRef<size_t>
     (temp_friends.first,
      location_pointer + 2,
      triangle));
  const facet facet2
    (TripleConstRef<size_t>
     (outer.second,
      outer.third,
      index),
     TripleConstRef<size_t>
     (temp_friends.second,
      triangle,
      location_pointer + 1));
  f.insert(f.end(), {facet1, facet2});
}

bool Delaunay::is_point_inside_big_triangle
(size_t index) const
{
  return InTriangle
    (TripleConstRef<Vector2D>
     (cor[olength],
      cor[olength + 1],
      cor[olength + 2]),
     cor[index]);
}

void Delaunay::add_point(size_t index,stack<std::pair<size_t, size_t> > &flip_stack)
{
  assert(is_point_inside_big_triangle(index));
  const size_t triangle = Walk(index);
  const Triplet<size_t> temp_friends(f[triangle].neighbors);
  update_f_in_add_point(triangle, temp_friends, index);
  update_friends_of_friends(triangle, temp_friends);

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
		{
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
	}

  check_if_flipping_is_needed(triangle, temp_friends, flip_stack);

  // _update number of facets
  location_pointer += 2;
}

void Delaunay::flip(size_t i, size_t j, stack<std::pair<size_t, size_t> > &flip_stack)
{
	if (j == last_loc)
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
		if (incircle(cor[prefetch_1.vertices.first],
			cor[prefetch_1.vertices.second],
			cor[prefetch_1.vertices.third],
			cor[check.first]) > 0)
		{
			//The point is in a circle change the facets and their friends
		  const auto v1 = prefetch_1.vertices[(other.second + 1) % 3];
			const auto f1 = prefetch_1.neighbors[other.second];
			const auto f12 = prefetch_1.neighbors[(other.second + 2) % 3];
			facet& prefetch_2 = f[indexes.second];
			const auto v2 = prefetch_2.vertices[(check.second + 1) % 3];
			const auto f2 = prefetch_2.neighbors[(check.second + 2) % 3];
			const auto f22 = prefetch_2.neighbors[check.second];
			prefetch_1.vertices.set(other.first, v1, check.first);
			prefetch_2.vertices.set(check.first, v2, other.first);
			prefetch_1.neighbors.set(f1, f2, indexes.second);
			prefetch_2.neighbors.set(f22, f12, indexes.first);
			// change the friends of the friends if needed
			if (f2 != last_loc)
			{
				f[f2].neighbors[find_index(f[f2], indexes.second)] = indexes.first;
			}
			if (f12 != last_loc)
			{
				f[f12].neighbors[find_index(f[f12], indexes.first)] = indexes.second;
			}
			// Calculate the new radius if needed
			if (CalcRadius)
			{
				radius[indexes.first] = CalculateRadius(indexes.first);
				radius[indexes.second] = CalculateRadius(indexes.second);
			}
			// clear the checked facets
			flip_stack.pop();
			// push into the stack the new facets to check
			if (prefetch_2.neighbors.first != last_loc)
			  flip_stack.push({indexes.second, prefetch_2.neighbors.first});
			if (prefetch_1.neighbors.second != last_loc)
				flip_stack.push(std::pair<size_t, size_t>(indexes.first,
					prefetch_1.neighbors.second));
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
	//	DataOnlyForBuild data;
	lastFacet = 0;
	CalcRadius = false;
	length = int(vp.size() + 3);
	size_t len = length - 3;
	olength = len;
	f.clear();
	cor.clear();
	f.reserve(2 * length + 1 + static_cast<size_t>(17 * sqrt(1.*length)));
	cor.reserve(length + 9 * static_cast<size_t>(sqrt(1.*length)));
	last_loc = INT_MAX;
	for (size_t i = 0; i < len; i++)
	{
		cor.push_back(vp[i]);
	}
	// Check point input
	CheckInput();

	// add the 3 extreme points
	const auto cellsize = CellSize(cell_points);
	const double width = cellsize[1] - cellsize[0];
	const double height = cellsize[3] - cellsize[2];
	const double extent = max(width, height);
	cor.insert(cor.end(),
		   {Vector2D(cellsize[0] - 100 * extent,
			     cellsize[2] - 100 * extent),
		       Vector2D(cellsize[1] + 100 * extent,
				cellsize[2] - 100 * extent),
		       Vector2D((cellsize[0] + cellsize[1]) / 2.0,
				cellsize[3] + 100 * extent)});
	// Create the big triangle, and assign friends
	f.push_back([&len, this]{
	    facet ans;
	    for(size_t i=0;i<3;++i){
	      ans.vertices[i] = len+i;
	      ans.neighbors[i] = last_loc;
	    }
	    return ans;}());
	location_pointer = 0;
	// add the points
	size_t nloop = length - 3;
	stack<std::pair<size_t, size_t> > flip_stack;
	for (size_t i = 0; i < nloop; i++)
		add_point(i,flip_stack);
	// Calculate radius
	radius.resize(f.size());
	for (size_t i = 0, n = f.size(); i < n; ++i)
		radius[i] = CalculateRadius(i);
	CalcRadius = true;
}

double Delaunay::triangle_area(size_t index)
{
	const TripleConstRef<Vector2D> p
		(cor[f[index].vertices.first],
		 cor[f[index].vertices.second],
		 cor[f[index].vertices.third]);
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

  boost::optional<size_t> walk_condition
  (const vector<Vector2D>& cor,
   const Triplet<size_t>& vertices,
   size_t point)
	{
	  for(size_t i=0;i<3;++i){
	    if (orient2d(TripleConstRef<Vector2D>
			 (cor[vertices[i]],
			  cor[vertices[(i+1)%3]],
			  cor[point])) < 0)
	      return i;
	  }
	  return boost::none;
	}

	size_t find_new_facet(const vector<Vector2D>& cor,
		const vector<facet>& f,
		size_t point,
		size_t last_facet)
	{
		size_t res = last_facet;
		auto next = walk_condition
		  (cor,
		   f[res].vertices,
		   point);
		while (next) {
		  res = f[res].neighbors[*next];
		  next = walk_condition(cor,
					f[res].vertices,
					point);
		}
		return res;
	}
}

size_t Delaunay::Walk(size_t point)
{
	lastFacet = find_new_facet(cor, f, point, lastFacet);
	return lastFacet;
}

vector<size_t> Delaunay::FindContainingTetras(size_t StartTetra, size_t point)
{
  vector<size_t> res;
  FindContainingTetras(StartTetra, point, res);
  return res;
}

double Delaunay::FindMaxRadius(size_t point)
{
	const vector<size_t> vec = FindContainingTetras(Walk(point), point);
	double r = 0;
	/*
	for (size_t i = 0; i < vec.size(); ++i)
		r = max(r, radius[vec[i]]);
	*/
	for(size_t vv : vec)
	  r = max(r, radius[vv]);
	return 2 * r;
}

void Delaunay::FindContainingTetras(size_t StartFacet, size_t point, vector<size_t> &result)
{
	result.clear();
	size_t PointLocation = FindPointInFacet(StartFacet, point);
	auto NextFacet = f[StartFacet].neighbors[PointLocation];
	result.reserve(12);
	result.push_back(NextFacet);
	while (NextFacet != StartFacet)
	{
	  PointLocation = FindPointInFacet(NextFacet, point);
		NextFacet = f[NextFacet].neighbors[PointLocation];
		result.push_back(NextFacet);
	}
}

size_t Delaunay::FindPointInFacet(size_t facet, size_t point)
{
	for (size_t i = 0; i < 3; ++i)
	  if (f[facet].vertices[i] == point)
			return i;
	UniversalError eo("Error in Delaunay, FindPointInFacet");
	eo.addEntry("Facet number", static_cast<int>(facet));
	eo.addEntry("Point number", static_cast<int>(point));
	throw eo;
}

bool Delaunay::IsOuterFacet(size_t facet)const
{
	//int PointNum=length-1;
	for (size_t i = 0; i<3; ++i)
		for (size_t j = 0; j < 3; ++j)
			if (f[facet].vertices[i] == olength + j)
				return true;
	return false;
}

double Delaunay::CalculateRadius(size_t facet)
{
	const double big = 1e10;
	std::array<double, 3> sides;
	for(int i=0;i<3;++i)
	  sides[i] = get_facet_coordinates(facet,i).distance
	    (get_facet_coordinates(facet,(i+1)%3));
	std::array<double, 3> temps;
	for(size_t i=0;i<3;++i){
	  temps[i] = sides[(i+1)%3]+sides[(i+2)%3]-sides[i];
	  if(temps[i] <= 0)
	    return (sides[i] > big*sides[(i+1)%3] ||
		    sides[i] > big*sides[(i+2)%3]) ?
	      0.5*sides[i] : 0.5*(sides[(i+1)%3]+sides[(i+2)%3]);
	}
	return sides[0]*sides[1]*sides[2]/
	  sqrt((sides[0]+sides[1]+sides[2])*
	       temps[0]*temps[1]*temps[2]);
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

double Delaunay::GetFacetRadius(size_t facet) const
{
	return radius[facet];
}

void Delaunay::ChangeOlength(size_t n)
{
  olength = n;
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

const facet& Delaunay::get_facet(size_t index) const
{
	return f[index];
}

const Vector2D& Delaunay::get_facet_coordinates
(size_t Facet, size_t vertice)
{
  return cor[f[Facet].vertices[vertice]];
}

const Vector2D& Delaunay::get_point(size_t index) const
{
	return cor[index];
}

size_t Delaunay::get_num_facet(void) const
{
	return f.size();
}

size_t Delaunay::get_length(void) const
{
	return length - 3;
}

size_t Delaunay::get_last_loc(void) const
{
  return last_loc;
}

void Delaunay::set_point(size_t index, Vector2D p)
{
  cor[index] = p;
}

size_t Delaunay::GetOriginalLength(void) const
{
  return olength;
}

vector<Vector2D>& Delaunay::GetMeshPoints(void)
{
	return cor;
}

size_t Delaunay::GetTotalLength(void)
{
	return cor.size();
}

void Delaunay::AddBoundaryPoints(vector<Vector2D> const& points)
{
	size_t n = points.size();
	stack<std::pair<size_t, size_t> > flip_stack;
	//	vector<int> order=HilbertOrder(points,n);
	for (size_t i = 0; i < n; ++i)
	{
		cor.push_back(points[i]);
		add_point(cor.size() - 1,flip_stack);
	}
}

void Delaunay::AddAditionalPoint(Vector2D const& vec)
{
	cor.push_back(vec);
}

size_t Delaunay::GetCorSize(void)const
{
	return cor.size();
}

bool Delaunay::IsTripleOut(size_t index) const
{
	int counter = 0;
	for(auto neighbor : f[index].neighbors)
	  counter += IsOuterFacet(neighbor);
	return counter > 1;
}

size_t Delaunay::FindTripleLoc(const facet& fct) const
{
	for (size_t i = 0; i < 3; ++i)
	  if (!IsOuterFacet(fct.neighbors[i]))
			return (i + 1) % 3;
	throw UniversalError("Trouble in constructing boundary triangles. No inner neighbor");
}

namespace
{
	bool IsOuterQuick(const facet& f, size_t olength)
	{
		for (size_t i = 0; i < 3; ++i)
		  if (f.vertices[i] >= olength)
				return true;
		return false;
	}

	bool IsEdgeFacet(const vector<facet>& facets, const facet& f, size_t olength)
	{
		for (size_t i = 0; i < 3; ++i)
		{
		  if (f.vertices[i] >= olength)
		    return false;
		  if (IsOuterQuick(facets[f.neighbors[i]], olength))
		    return true;
		}
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

vector<size_t> Delaunay::GetOuterFacets(size_t start_facet, size_t real_point, size_t olength2)
{
	size_t cur_facet = start_facet;
	vector<size_t> f_temp, containing_facets;
	f_temp.reserve(static_cast<size_t>(10 * sqrt(static_cast<double>(olength2))));
	size_t point_index = FindPointInFacet(cur_facet, real_point);
	if (IsOuterQuick(f[f[cur_facet].neighbors[point_index]], olength2))
	{
		point_index = (point_index + 1) % 3;
		real_point = f[cur_facet].vertices[point_index];
	}
	do
	{
		FindContainingTetras(cur_facet, real_point, containing_facets);
		size_t old_current = cur_facet;
		for (size_t i = 0; i < containing_facets.size(); ++i)
		{
		  if (IsEdgeFacet(f, f[(containing_facets[i])], olength2) &&
		      containing_facets[i] != old_current)
		    cur_facet = containing_facets[i];
			if (!IsOuterQuick(f[containing_facets[i]], olength2))
			  f_temp.push_back(containing_facets[i]);
		}
		point_index = (1 + FindPointInFacet(cur_facet, real_point)) % 3;
		if (IsTripleOut(cur_facet))
			point_index = (point_index + 1) % 3;
		real_point = f[cur_facet].vertices[point_index];
	} while (start_facet != cur_facet);
	sort(f_temp.begin(), f_temp.end());
	f_temp = unique(f_temp);
	return f_temp; //adapter1<size_t,int>(f_temp);
}

vector<vector<size_t> > Delaunay::FindOuterPoints(vector<Edge> const& edges)
{
	// We add the points in a counter clockwise fashion
	if (olength < 100)
	  return vector<vector<size_t> >
	    (edges.size(),
	     [](size_t n)
	     {vector<size_t> ans(n);
	       iota(ans.begin(), ans.end(),0);
	       return ans;}(olength));
	vector<vector<int> > res(edges.size());
	vector<int> res_temp, outer_points, f_temp;
	res_temp.reserve(static_cast<size_t>(20 * sqrt(1.0*static_cast<double>(olength))));
	f_temp.reserve(static_cast<size_t>(10 * sqrt(1.0*static_cast<double>(olength))));
	outer_points.reserve(static_cast<size_t>(10 * sqrt(1.0*static_cast<double>(olength))));
	// Walk to an outer point
	size_t cur_facet = Walk(olength);
	vector<vector<size_t> > toduplicate(edges.size());
	vector<bool> checked(f.size(), false);
	AddOuterFacets(cur_facet, toduplicate, edges, checked);
	for (size_t i = 0; i < edges.size(); ++i)
	{
		sort(toduplicate[i].begin(), toduplicate[i].end());
		toduplicate[i] = unique(toduplicate[i]);
	}
	return toduplicate;
}

void Delaunay::AddRigid(vector<Edge> const& edges,
	vector<vector<size_t> > &toduplicate)
{
	vector<size_t> toremove;
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
			Vector2D temp = cor[toduplicate[i][j]] - edges[i].vertices.first;
			temp = 2 * par*ScalarProd(par, temp) - temp + edges[i].vertices.first;
			if (InTriangle(TripleConstRef<Vector2D>
				(cor[olength],
					cor[olength + 1],
					cor[olength + 2]),
				temp))
				toadd.push_back(temp);
			else
				toremove.push_back(j);
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
			eo.addEntry("Error in AddRigid", 0);
			throw;
		}
		ReArrangeVector(toduplicate[i], order);
	}
}

namespace
{
	vector<Edge> GetCornerEdges(OuterBoundary const& obc)
	{
		const double dx = obc.GetGridBoundary(Right) - obc.GetGridBoundary(Left);
		const double dy = obc.GetGridBoundary(Up) - obc.GetGridBoundary(Down);
		vector<Edge> res;
		const Vector2D RU(obc.GetGridBoundary(Right), obc.GetGridBoundary(Up));
		const Vector2D LU(obc.GetGridBoundary(Left), obc.GetGridBoundary(Up));
		const Vector2D LD(obc.GetGridBoundary(Left), obc.GetGridBoundary(Down));
		const Vector2D RD(obc.GetGridBoundary(Right), obc.GetGridBoundary(Down));
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

vector<vector<size_t> > Delaunay::AddPeriodic(const OuterBoundary& obc, vector<Edge> const& edges,
	vector<vector<size_t> > &toduplicate)
{
	const double dx = obc.GetGridBoundary(Right) - obc.GetGridBoundary(Left);
	const double dy = obc.GetGridBoundary(Up) - obc.GetGridBoundary(Down);
	const std::array<Vector2D,4> changes1 = 
	  {Vector2D(-dx,0),
	   Vector2D(0,-dy),
	   Vector2D(dx,0),
	   Vector2D(0,dy)};
	for (size_t i = 0; i < edges.size(); ++i)
	{
		if (toduplicate[i].empty())
			continue;
		const Vector2D change = changes1[i];
		vector<Vector2D> toadd(toduplicate[i].size());
		transform(toduplicate[i].begin(),
			  toduplicate[i].end(),
			  toadd.begin(),
			  [&change,this](size_t x)
			  {return cor[x]+change;});
		vector<int> order = HilbertOrder(toadd, static_cast<int>(toadd.size()));
		ReArrangeVector(toadd, order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(toduplicate[i], order);
		//toduplicate[i]=pointstemp;
	}
	// Done with sides do corners now
	vector<Edge> corneredges = GetCornerEdges(obc);
	vector<vector<size_t> > corners(toduplicate.size());
	for (size_t i = 0; i < toduplicate.size(); ++i){
	  for (size_t j = 0; j < toduplicate[i].size(); ++j){
	      const size_t facet_loc = Walk(toduplicate[i][j]);
	      const Vector2D center = cor[toduplicate[i][j]];
	      const double R = 2*GetMaxRadius(toduplicate[i][j], facet_loc);
	      if (CircleSegmentIntersect(corneredges[2*i], center, R))
		corners[i].push_back(toduplicate[i][j]);
	      if (CircleSegmentIntersect(corneredges[(2*i+7)%8], center, R))
		corners[(i + 3) % 4].push_back(toduplicate[i][j]);
	    }
		 }
	  const std::array<Vector2D, 4> changes = 
	    {Vector2D(-dx,-dy),
	     Vector2D(dx,-dy),
	     Vector2D(dx,dy),
	     Vector2D(-dx,dy)};
	for (size_t i = 0; i < corners.size(); ++i)
	{
		if (corners[i].empty())
			continue;
		sort(corners[i].begin(), corners[i].end());
		corners[i] = unique(corners[i]);
		const Vector2D change = changes[i];
		vector<Vector2D> toadd;
		toadd.reserve(corners[i].size());
		for(auto x : corners[i])
		  toadd.push_back(cor[x]+change);
		vector<int> order = HilbertOrder(toadd, static_cast<int>(toadd.size()));
		ReArrangeVector(toadd, order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(corners[i], order);
	}
	return corners;
}

void Delaunay::AddHalfPeriodic(OuterBoundary const& obc, vector<Edge> const& edges,
	vector<vector<size_t> > &toduplicate)
{
	const double dx = obc.GetGridBoundary(Right) - obc.GetGridBoundary(Left);
	//	const double dy=obc.GetGridBoundary(Up)-obc.GetGridBoundary(Down);
	const std::array<Vector2D,4> changes = 
	  {Vector2D(-dx,0),
	   Vector2D(-dx,0),
	   Vector2D(dx,0),
	   Vector2D(dx,0)};
	for (size_t i = 0; i < edges.size(); ++i)
	{
		if (toduplicate[i].empty())
			continue;
		const Vector2D change = changes[i];
		vector<Vector2D> toadd;
		toadd.reserve(toduplicate[i].size());
		//vector<int> pointstemp(toduplicate[i].size());
		Vector2D par(Parallel(edges[i]));
		par = par / abs(par);
		for (size_t j = 0; j < toduplicate[i].size(); ++j)
		{
			Vector2D temp = cor[toduplicate[i][j]];
			if (i % 2 == 1)
			{
				temp -= edges[i].vertices.first;
				temp = 2 * par*ScalarProd(par, temp) - temp + edges[i].vertices.first;
			}
			toadd.push_back(temp + change);
			//pointstemp[j]=j;
		}
		       vector<size_t> order = adapter1<int,size_t>(HilbertOrder(toadd, static_cast<int>(toadd.size())));
		ReArrangeVector(toadd, order);
		AddBoundaryPoints(toadd);
		ReArrangeVector(toduplicate[i], order);
		//toduplicate[i]=pointstemp;
	}
}

vector<vector<size_t> > Delaunay::BuildBoundary(const OuterBoundary& obc, vector<Edge> const& edges)
{
  vector<vector<size_t> > toduplicate = FindOuterPoints(edges);
#ifdef RICH_MPI
	OrgIndex.clear();
#endif
	if (obc.GetBoundaryType() == Rectengular)
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
		if (obc.GetBoundaryType() == Periodic)
		{
			vector<vector<size_t> > corners = AddPeriodic(obc, edges, toduplicate);
			for (size_t i = 0; i < 4; ++i)
				toduplicate.push_back(corners[i]);
		}
		else
		{
			AddHalfPeriodic(obc, edges, toduplicate);
		}
	}
	return toduplicate;
}

Vector2D Delaunay::GetCircleCenter(size_t index)const
{
	const facet& F = f[index];
	double x1 = cor[F.vertices[0]].x;
	double x2 = cor[F.vertices[1]].x;
	double x3 = cor[F.vertices[2]].x;
	double y1 = cor[F.vertices[0]].y;
	double y2 = cor[F.vertices[1]].y;
	double y3 = cor[F.vertices[2]].y;
	// Do we have a case where two point are very close compared to the third?
	const double d12 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
	const double d23 = (x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2);
	const double d13 = (x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3);
	if(d23 < 0.1*(d13+d12)){
	  x1 -= x2;
	  x3 -= x2;
	  y1 -= y2;
	  y3 -= y2;
	  const double d_inv = 1 / (2 * (x3*y1 - y3*x1));
	  return Vector2D
	    ((y1*(x3*x3 + y3*y3) - y3*(x1*x1 + y1*y1))*d_inv + x2,
	     (x3*(x1*x1 + y1*y1) - x1*(x3*x3 + y3*y3))*d_inv + y2);
	}
	x2 -= x1;
	x3 -= x1;
	y2 -= y1;
	y3 -= y1;
	const double d_inv = 1 / (2 * (x2*y3 - y2*x3));
	return Vector2D
	  ((y3*(x2*x2 + y2*y2) - y2*(x3*x3 + y3*y3))*d_inv + x1,
	   (-x3*(x2*x2 + y2*y2) + x2*(x3*x3 + y3*y3))*d_inv + y1);
}

double Delaunay::GetMaxRadius(size_t point, size_t startfacet)
{
	double res = 0;
	const vector<size_t> neigh = FindContainingTetras(startfacet, point);
	for (size_t ni : neigh)
	  res = max(res, radius[ni]);
	return res;
}

void Delaunay::AddOuterFacets(size_t tri, vector<vector<size_t> > &toduplicate,
	vector<Edge> const& edges, vector<bool> &checked)
{
	stack<size_t> tocheck;
	tocheck.push(tri);
	while (!tocheck.empty())
	{
		size_t cur_facet = tocheck.top();
		tocheck.pop();
		for (size_t i = 0; i < 3; ++i)
		{
			bool added = false;
			if (checked[f[cur_facet].vertices[i]] || 
			    f[cur_facet].vertices[i] >= olength)
				continue;
			vector<size_t> neigh = FindContainingTetras(cur_facet, f[cur_facet].vertices[i]);
			for (size_t k = 0; k < neigh.size(); ++k)
			{
				Vector2D center = GetCircleCenter(neigh[k]);
				for (size_t l = 0; l < edges.size(); ++l)
				{
					if (CircleSegmentIntersect(edges[l], center, radius[neigh[k]]))
					{
					  toduplicate[l].push_back(f[cur_facet].vertices[i]);
						added = true;
					}
				}
			}
			checked[f[cur_facet].vertices[i]] = true;
			if (added)
			  for_each(neigh.begin(),
				   neigh.end(),
				   [&tocheck](size_t x)
				   {tocheck.push(x);});
		}
	}
}

#ifdef RICH_MPI
size_t Delaunay::findSomeOuterPoint(void)
{
	const size_t cur_facet = Walk(olength);
	for (size_t i = 0; i < 3; ++i) {
	  const size_t candidate = f.at(cur_facet).vertices[i];
	  if (candidate < olength)
	      return candidate;
	}
	assert(false && "something went wrong");
	return 0;
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

	stack<size_t> initialise_tocheck
		(const vector<size_t>& neightemp)
	{
		stack<size_t> res;
		for (size_t i = 0; i < neightemp.size(); ++i)
			res.push(neightemp[i]);
		return res;
	}

	vector<size_t> calc_self_intersection
	(const vector<Edge>& edge_list,
	 const Circle& circle)
	{
		vector<size_t> res;
		for (size_t i = 0; i < edge_list.size(); ++i) {
			const Edge& edge = edge_list.at(i);
			if (edge_circle_intersect(edge, circle))
				res.push_back(i);
		}
		return res;
	}
}

vector<vector<size_t> > Delaunay::AddOuterFacetsMPI
(int point,
	vector<vector<size_t> > &toduplicate,
	vector<int> &neigh,
	vector<bool> &checked,
	Tessellation const &tproc,
	const vector<Edge>& own_edges,
	bool recursive)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<vector<size_t> > res;
	vector<int> vtemp;
	if (!recursive)
		res.resize(own_edges.size());
	stack<size_t> tocheck = initialise_tocheck
	  (FindContainingTetras(Walk(point), point));
	if (recursive)
	{
		vector<size_t> allouter;
		for (size_t i = 0; i < toduplicate.size(); ++i)
		{
			for (size_t j = 0; j < toduplicate[i].size(); ++j)
			{
			  vector<size_t> temp = FindContainingTetras(Walk(toduplicate[i][j]), toduplicate[i][j]);
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
		size_t cur_facet = tocheck.top();
		tocheck.pop();
		for (size_t i = 0; i < 3; ++i)
		{
			bool added = false;
			size_t max_neigh = 0;
			if (f[cur_facet].vertices[i] >= olength ||
			    checked[f[cur_facet].vertices[i]])
				continue;
			vector<int> neighs = adapter1<size_t,int>(FindContainingTetras(cur_facet, f[cur_facet].vertices[i]));
			for (size_t k = 0; k < neighs.size(); ++k)
			{
				Circle circ(GetCircleCenter(neighs[k]), radius[neighs[k]]);
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
				  const vector<size_t> self_intersection = calc_self_intersection(own_edges, circ);
					if (!self_intersection.empty())
						added = true;
					BOOST_FOREACH(size_t sindex, self_intersection)
					  res[sindex].push_back(f[cur_facet].vertices[i]);
				}
				else
				{
					for (size_t jj = 0; jj < 3; ++jj)
					  max_neigh = max(max_neigh, f[neighs[k]].vertices[jj]);
				}
				if (!cputosendto.empty())
				{
					added = true;
					for (size_t j = 0; j < cputosendto.size(); ++j)
					{
					  size_t index = find(neigh.begin(), neigh.end(), cputosendto[j]) - neigh.begin();
						if (index < neigh.size())
						  toduplicate.at(index).push_back(f[cur_facet].vertices[i]);
						else {
							neigh.push_back(cputosendto[j]);
							toduplicate.push_back
							  (vector<size_t>(1, f[cur_facet].vertices[i]));
						}
					}
				}
			}
			checked[f[cur_facet].vertices[i]] = true;
			if (added||(recursive && max_neigh>=olength))
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

pair<vector<vector<size_t> >, vector<vector<size_t> > >
Delaunay::findOuterPoints
(const Tessellation& t_proc,
	const vector<Edge>& edge_list,
	const vector<Edge>& box_edges,
	vector<vector<int> > &NghostIndex)
{
	auto neighbors_own_edges =
		calc_neighbors_own_edges(t_proc, edge_list);
	const size_t some_outer_point = findSomeOuterPoint();

	vector<vector<size_t> > to_duplicate(neighbors_own_edges.size());
	vector<bool> checked(olength, false);
	vector<vector<size_t> > self_points =
		AddOuterFacetsMPI
		(static_cast<int>(some_outer_point),
			to_duplicate, // indices of points to send
		 neighbors_own_edges, // Rank of processes to send to
			checked,
			t_proc,
			box_edges);
	BOOST_FOREACH(vector<size_t>& line, to_duplicate) 
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
			size_t location = std::find(neighbors_own_edges.begin(), neighbors_own_edges.end(),
				status.MPI_SOURCE) - neighbors_own_edges.begin();
			if (location >= neighbors_own_edges.size())
				throw UniversalError("Bad location in mpi exchange");
			try
			{
				incoming[location] = list_unserialize(temprecv, cor[0]);
			}
			catch (UniversalError &eo)
			{
				eo.addEntry("Error in first send in triangulation", 0.0);
				throw;
			}
		}
		else
			if (status.MPI_TAG != 1)
				throw UniversalError("Wrong mpi tag");
	}
	MPI_Waitall(static_cast<int>(req.size()), &req[0], MPI_STATUSES_IGNORE);
	// Incorporate points recieved into triangulation
	BOOST_FOREACH(vector<size_t> &line, self_points)
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
	return pair<vector<vector<size_t> >, vector<vector<size_t> > >
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

vector<vector<size_t> >
Delaunay::boundary_intersection_check
(const vector<Edge>& edges,
	const vector<vector<size_t> >& to_duplicate)
{
	vector<vector<size_t> > res;
	res.reserve(edges.size());
	BOOST_FOREACH(const Edge& edge, edges) {
		res.push_back(vector<size_t>());
		BOOST_FOREACH(const vector<size_t>& line, to_duplicate) {
			BOOST_FOREACH(const size_t index, line) {
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
	vector<vector<size_t> > &to_duplicate,
	vector<vector<int> >& self_points,
	const vector<Edge>& box_edges,
	vector<vector<size_t> >& NghostIndex)
{
	const vector<vector<size_t> > boundary_points =
		boundary_intersection_check(box_edges, to_duplicate);
	vector<vector<size_t> > real_boundary_points
		(box_edges.size());
	for (size_t i = 0; i < boundary_points.size(); ++i)
	{
		sort(self_points.at(i).begin(), self_points.at(i).end());
		BOOST_FOREACH(size_t bp, boundary_points.at(i))
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

	vector<vector<size_t> > to_duplicate_2 = to_duplicate;
	BOOST_FOREACH(vector<size_t>& line, to_duplicate_2)
		sort(line.begin(), line.end());
	vector<int> neighbors_own_edges =
		calc_neighbors_own_edges(t_proc, edge_list);
	vector<int> old_neighbors = neighbors_own_edges;
	assert(!to_duplicate.empty());
	BOOST_FOREACH(const vector<size_t>& line, to_duplicate)
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
				const size_t my_index = it - old_neighbors.begin();
				if (!binary_search
					(to_duplicate_2.at(my_index).begin(),
						to_duplicate_2.at(my_index).end(),
						to_duplicate.at(i).at(j)))
				  messages.at(i).push_back(static_cast<int>(to_duplicate.at(i).at(j)));
			}
			else
			  messages.at(i).push_back(static_cast<int>(to_duplicate.at(i).at(j)));
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
		auto size = static_cast<int>(tosend[i].size());
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
			size_t location = std::find(neighbors_own_edges.begin(), neighbors_own_edges.end(), status.MPI_SOURCE) - neighbors_own_edges.begin();
			if (location >= neighbors_own_edges.size())
				throw UniversalError("Bad location in mpi exchange");
			try
			{
				incoming[location] = list_unserialize(temprecv, cor[0]);
			}
			catch (UniversalError &eo)
			{
				eo.addEntry("Error in second send in triangulation", 0.0);
				eo.addEntry("Mpi status", static_cast<double>(status.MPI_SOURCE));
				eo.addEntry("Mpi tag", static_cast<double>(status.MPI_TAG));
				eo.addEntry("Mpi count", static_cast<double>(count));
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
			NghostIndex[i].push_back(cor.size() + j);
			OrgIndex.push_back(cor.size() + j);
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

	return  pair<vector<vector<int> >,vector<int> > 
	  (adapter2<size_t,int>(to_duplicate),
	   neighbors_own_edges);
	//	   adapter1<size_t,int>(neighbors_own_edges));
}

pair<vector<vector<int> >, vector<int> > Delaunay::BuildBoundary
(const OuterBoundary& obc,
 const Tessellation& tproc,
 vector<vector<int> >& Nghost)
{
	vector<Edge> edges;
	OrgIndex.clear();
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<int> edge_index = tproc.GetCellEdges(rank);
	for (size_t i = 0; i < edge_index.size(); ++i)
		edges.push_back(tproc.GetEdge(edge_index[i]));
	vector<Edge> box_edges = obc.GetBoxEdges();
	auto to_duplicate =
		findOuterPoints(tproc, edges, box_edges, Nghost);
	auto aux1 = adapter2<int, size_t>(Nghost);
	auto aux2 = adapter2<size_t, int>(to_duplicate.second);
	auto ans = FindOuterPoints2(tproc,edges,
				    to_duplicate.first,
				    aux2,
				    box_edges, aux1);
	Nghost = adapter2<size_t, int>(aux1);
	to_duplicate.second = adapter2<int, size_t>(aux2);
	return ans;
}

size_t Delaunay::GetOrgIndex(size_t index)const
{
	if (index < olength)
		return olength;
	else
		return OrgIndex.at(index - 3 - olength);
}

#endif // RICH_MPI
