// TryVoro++.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <string>
#include <voro++.hh>

using namespace std;

int main()
{
	cout << "Hello, TryVoro++" << endl;

	// Create a 4x4x4 set of points, starting at with indices in (-1.5, -0.5, +0.5, +1.5)
	voro::container con(-2, 2,  // x
		-2, 2, // y
		-2, 2, // z
		6, 6, 6, // Blocks
		false, false, false, 8);

	int n = 0;
	for (double x = -1.5; x <= 1.5; x += 1.0)
		for (double y = -1.5; y <= 1.5; y += 1.0)
			for (double z = -1.5; z <= 1.5; z += 1.0)
			{
				con.put(n, x, y, z);
			}

	voro::c_loop_all looper(con);
	for (looper.start(); looper.inc();)
	{
		voro::voronoicell cell;
		con.compute_cell(cell, looper);
		cell.output_face_vertices();
		cout << endl;
	}

	string s;
	cin >> s;
}