#include <iostream>
#include <fstream>
#include <cmath>
#include "source/newtonian/two_dimensional/spatial_distributions/triangle.hpp"

using namespace std;

int main(void)
{
  vector<Vector2D> vv(3);
  vv[0].Set(0.5,0.6);
  vv[1].Set(0.7,0.5);
  vv[2].Set(0.4,0.4);
  Triangle t(vv,1,0);
  
  ofstream f;
  Vector2D p1(0.5,0.5);
  Vector2D p2(1,1);
  f.open("res.txt");
  f << t.EvalAt(p1) << endl;
  f << t.EvalAt(p2) << endl;
  f.close();
}
