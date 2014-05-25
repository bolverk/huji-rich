#include <iostream>
#include <fstream>
#include "source/newtonian/common/hydrodynamic_variables.hpp"

// This test checks that the subscripting in the primitive class is working properly

using namespace std;

int main(void)
{
  Primitive p;
  p[0] = 1.0;
  p.Pressure = 2.0;
  p.Energy = 3.0;
  p.SoundSpeed = 4.0;
  p.Velocity.x = 5.0;
  p.Velocity.y = 6.0;

  ofstream f("res.txt");
  for(int i=0;i<6;i++)
    f << p[i] << endl;
  f.close();
}
