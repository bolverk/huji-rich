#include <iostream>
#include <fstream>
#include "source/newtonian/common/hllc.hpp"

// This test checks that the solver gives the right result for hydrostatic conditions

using namespace std;

int main(void)
{
  Primitive left;
  Primitive right;
  right.Density = left.Density = 1.2;
  right.Pressure = left.Pressure = 3.4;
  right.Energy = left.Energy = 5.6;
  right.SoundSpeed = left.SoundSpeed = 7.8;
  right.Velocity.x = left.Velocity.x = 0;
  right.Velocity.y = left.Velocity.y = 0;
  double velocity = 9.1;
  Hllc solver;
  Conserved res = solver(left, right, velocity);

  ofstream f("res.txt");
  f << res.Mass << endl;
  f << -right.Density*velocity << endl;
  f.close();
}
