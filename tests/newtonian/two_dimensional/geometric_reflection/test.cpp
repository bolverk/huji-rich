#include <iostream>
#include <fstream>
#include "source/tessellation/geometry.hpp"

using namespace std;

int main(void)
{
  Vector2D v(1,1);
  Vector2D axis(0,2);
  Vector2D res = Reflect(v, axis);
  ofstream f("res.txt");
  f << res.x << " "
    << res.y << endl;
  f.close();
}
