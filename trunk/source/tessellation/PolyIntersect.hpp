#include "geotests.hpp"
#include "Edge.hpp"

enum InFlags {UnKnown,Pi,Qi};
enum IntersectFlags {True,False,Par};

vector<Vector2D> ConvexIntersect(vector<Vector2D> const& poly0,vector<Vector2D>
	const& poly1);

IntersectFlags SegmentIntersection(Vector2D const& p0,Vector2D const& p1,
	Vector2D const& q0,Vector2D const& q1,Vector2D &Intersection);

vector<Vector2D> GetParEdge(Vector2D const& p0,Vector2D const& p1,
	Vector2D const& q0,Vector2D const& q1);