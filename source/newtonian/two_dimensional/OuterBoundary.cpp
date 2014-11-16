#include "OuterBoundary.hpp"

OuterBoundary::~OuterBoundary(void) {}

vector<Edge> OuterBoundary::GetBoxEdges(void) const
{
	vector<Edge> res(4);
	res[1]=Edge(Vector2D(GetGridBoundary(Left),GetGridBoundary(Up)),
		Vector2D(GetGridBoundary(Right),GetGridBoundary(Up)),0,0);
	res[3]=Edge(Vector2D(GetGridBoundary(Left),GetGridBoundary(Down)),
		Vector2D(GetGridBoundary(Right),GetGridBoundary(Down)),0,0);
	res[0]=Edge(Vector2D(GetGridBoundary(Right),GetGridBoundary(Up)),
		Vector2D(GetGridBoundary(Right),GetGridBoundary(Down)),0,0);
	res[2]=Edge(Vector2D(GetGridBoundary(Left),GetGridBoundary(Up)),
		Vector2D(GetGridBoundary(Left),GetGridBoundary(Down)),0,0);
	return res;
}

