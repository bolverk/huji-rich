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
	voro::container con(-1, 1,  // x
		-1, 1, // y
		-1, 1, // z
		6, 6, 6, // Blocks
		false, false, false, 8);

	int n = 0;
	for (double x = -0.5; x <= 0.5; x += 1.0)
		for (double y = -0.5; y <= 0.5; y += 1.0)
			for (double z = -0.5; z <= 0.5; z += 1.0)
			{
				con.put(n, x, y, z);
			}

	con.draw_cells_pov("cube.pov");

	voro::c_loop_all looper(con);
	for (looper.start(); looper.inc();)
	{
		voro::voronoicell cell;
		double x, y, z;
		looper.pos(x, y, z);

		cout << "Cell at (" << x << ", " << y << ", " << z << ")" << endl;

		con.compute_cell(cell, looper);
		vector<double> vertices;
		cell.vertices(vertices);

		for (int i = 0; i < vertices.size(); i += 3)
		{
			vertices[i] += x;
			vertices[i + 1] += y;
			vertices[i + 2] += z;
		}

		vector<int> faces, orders;
		cell.face_vertices(faces);
		cell.face_orders(orders);
		for (int i = 0; i < faces.size(); i++)
		{
			int v = faces[i] * 3;
			cout << "(" << vertices[v] << ", " << vertices[v + 1] << +", " << vertices[v + 2] << ") ";
		}
		cout << endl;
		cell.output_face_vertices();
		cout << endl;
	}

	string s;
	cin >> s;
}