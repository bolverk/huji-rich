#include <iostream>
#include "source/tessellation/geometry.hpp"
#include "source/tessellation/Delaunay.hpp"
#include "source/misc/simple_io.hpp"

using namespace std;

int main(void)
{
#ifdef RICH_MPI
  Delaunay tri;
#else

  write_number(0,"serial_ignore.txt");

#endif // RICH_MPI
  return 0;
}
