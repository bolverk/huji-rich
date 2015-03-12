// StressVoronoi.cpp : Defines the entry point for the console application.
//

#include "GeometryCommon/Vector3D.hpp"
#include "GeometryCommon/OuterBoundary3D.hpp"
#include "Voronoi/VoroPlusPlus.hpp"
#include "Voronoi/GhostBusters.hpp"
#include "Voronoi/TetGenDelaunay.hpp"
#include "Voronoi/TetGenTessellation.hpp"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include "VectorNamer.hpp"
#include "CommandLineParser.hpp"

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

namespace fs = boost::filesystem;
using namespace std;

static double RandomDouble(double min, double max);
static Vector3D RandomPoint(const OuterBoundary3D &boundary);
static void RandomBoundary(Vector3D &frontTopRight, Vector3D &backLowerLeft);
static vector<Vector3D> RandomPoints(int num, const OuterBoundary3D &boundary);

static bool Initialize();
static vector<Vector3D> ReadPoints(const std::string filename);
static void WritePoints(const std::vector<Vector3D> &points, const std::string &filename);
static void WritePoints(const std::vector<VectorRef> &points, const std::string &filename);
static void WriteMatlabPoints(const std::vector<VectorRef> &points, const std::string &filename);
static void WriteMatlabPoints(const std::vector<Vector3D> &points, const std::string &filename);

static void RunVoronoi(Tessellation3D *tes, std::string name);

static void CompareTessellations();
static void CompareTessellations(const std::string name1, const std::string name2);
static void CompareCells(const Tessellation3D *tes1, const Tessellation3D *tes2, size_t cellNum);
static bool IsSameFace(const Face &face1, const Face &face2);
static double FaceSimilarity(const Face &face1, const Face &face2);

static VectorNamer namer;
static Arguments args;
static OuterBoundary3D *boundary;
static vector<Vector3D> points;

static map<string, Tessellation3D *> tessellations;

int main(int argc, char*argv[])
{
	if (!ParseArguments(argc, argv, args))
		return -1;

	if (!Initialize())
		return -1;

	if (args.RunBruteForce)
	{
		TetGenTessellation<BruteForceGhostBuster> *del = new TetGenTessellation<BruteForceGhostBuster>();
		RunVoronoi(del, "brute-force");
		WritePoints(del->AllPoints, "brute-force.node");
		WriteMatlabPoints(del->AllPoints, "brute-force.m");
	}

	if (args.RunFullBruteForce)
	{
		TetGenTessellation<FullBruteForceGhostBuster> *del = new TetGenTessellation<FullBruteForceGhostBuster>();
		RunVoronoi(del, "full-brute-force");
		WritePoints(del->AllPoints, "full-brute-force.node");
		WriteMatlabPoints(del->AllPoints, "full-brute-force.m");
	}

	if (args.RunCloseToBoundary)
	{
		TetGenTessellation<CloseToBoundaryGhostBuster> *del = new TetGenTessellation<CloseToBoundaryGhostBuster>();
		RunVoronoi(del, "close-to-boundary");
		WritePoints(del->AllPoints, "close-to-boundary.node");
		WriteMatlabPoints(del->AllPoints, "close-to-boundary.m");
	}

	if (args.RunVoroPlusPlus)
	{
		VoroPlusPlus *voro = new VoroPlusPlus();
		RunVoronoi(voro, "voro++");
	}

	CompareTessellations();

	cout << endl << "Stress test done." << endl;

	for (map<string, Tessellation3D *>::iterator it = tessellations.begin(); it != tessellations.end(); it++)
		delete it->second;

#ifdef _DEBUG
	string s;
	getline(cin, s);
#endif
}

static bool Initialize()
{
	cout << "Initializing..." << endl;
	// srand(17);
	srand((unsigned int)time(NULL));

	namer.SaveZero();
	Vector3D ftr, bll;
	RandomBoundary(ftr, bll);
	boundary = new RectangularBoundary3D(ftr, bll);

	if (!args.InputFile.empty())
	{
		points = ReadPoints(args.InputFile);
		if (points.empty())
		{
			cerr << "Invalid node file " << args.InputFile << endl;
			return false;
		}
		cout << "Read " << points.size() << " points" << endl;
	}
	else
	{
		points = RandomPoints(args.NumPoints, *boundary);
		cout << "Generated " << points.size() << " points" << endl;
	}
	WritePoints(points, "orig.node");
	WriteMatlabPoints(points, "orig.m");

	for (vector<Vector3D>::iterator it = points.begin(); it != points.end(); it++)
		namer.GetName(*it, "C");
	return true;
}

static void WritePoints(const vector<VectorRef> &pointRefs, const std::string &filename)
{
	vector<Vector3D> points;
	points.reserve(pointRefs.size());

	for (vector<VectorRef>::const_iterator it = pointRefs.begin(); it != pointRefs.end(); it++)
		points.push_back(**it);

	WritePoints(points, filename);
}

static void WritePoints(const vector<Vector3D> &points, const std::string &filename)
{
	ofstream output;
	Tetrahedron big = Delaunay::CalcBigTetrahedron(*boundary);
	fs::path full_path = fs::path(args.OutputDirectory) / filename;
	output.open(full_path.string());
	output.precision(10);
	output << "# Nodes, Dim, #Attrs, Boundary" << endl;
	output << points.size() + 4 << " 3 0 1" << endl;
	output << "# Index, X, Y, Z, Boundary" << endl;
	output << "# Main points" << endl;

	for (size_t i = 0; i < points.size(); i++)
		output << i + 1 << " " << points[i].x << " " << points[i].y << " " << points[i].z << " 0" << endl;

	output << "# Boundary points" << endl;
	for (int i = 0; i < 4; i++)
		output << i + 1 + points.size() << " " << big[i]->x << " " << big[i]->y << " " << big[i]->z << " 1" << endl;

	output << endl;
}

static void WriteMatlabPoints(const vector<VectorRef> &pointRefs, const std::string &filename)
{
	vector<Vector3D> points;
	points.reserve(pointRefs.size());

	for (vector<VectorRef>::const_iterator it = pointRefs.begin(); it != pointRefs.end(); it++)
		points.push_back(**it);

	WriteMatlabPoints(points, filename);
}

static void WriteMatlabPoints(const vector<Vector3D> &points, const std::string &filename)
{
	ofstream output;
	Tetrahedron big = Delaunay::CalcBigTetrahedron(*boundary);
	fs::path full_path = fs::path(args.OutputDirectory) / filename;
	output.open(full_path.string());
	output.precision(10);
	
	output << "m = [";
	for (size_t i = 0; i < points.size(); i++)
	{
		if (i > 0)
			output << "," << endl;
		output << "[" << points[i].x << ", " << points[i].y << ", " << points[i].z << "]";
	}

	for (size_t i = 0; i < 4; i++)
	{
		output << "," << endl;
		output << "[" << big[i]->x << ", " << big[i]->y << ", " << big[i]->z << "]";
	}
	output << "]" << endl;
}

static vector<Vector3D> ReadPoints(const std::string filename)
{
	ifstream input;
	input.open(filename);

	vector<Vector3D> points;
	string line;

	bool first = true;
	while (getline(input, line))
	{
		boost::trim(line);
		if (line.empty())
			continue;
		if (line[0] == '#')
			continue;
		if (first)  // Ignore the first line, which is just the count
		{
			first = false;
			continue;
		}

		int index, boundary = 0;
		double x, y, z;
		stringstream strm(line);
		strm >> index >> x >> y >> z >> boundary;

		if (!boundary)
			points.push_back(Vector3D(x, y, z));
	}

	return points;
}

static double RandomDouble(double min, double max)
{
	double fraction = (double)rand() / RAND_MAX;
	return (max - min) * fraction + min;
}

static Vector3D RandomPoint(const OuterBoundary3D &boundary)
{
	Vector3D v;

	v.x = RandomDouble(boundary.BackLowerLeft().x, boundary.FrontUpperRight().x);
	v.y = RandomDouble(boundary.BackLowerLeft().y, boundary.FrontUpperRight().y);
	v.z = RandomDouble(boundary.BackLowerLeft().z, boundary.FrontUpperRight().z);

	return v;
}

static void RandomBoundary(Vector3D &frontTopRight, Vector3D &backLowerLeft)
{
	backLowerLeft.x = -1000;
	frontTopRight.x = -500;

	backLowerLeft.y = 0;
	frontTopRight.y = 1000;
	
	backLowerLeft.z = -300;
	frontTopRight.z = 300;
}

static vector<Vector3D> RandomPoints(int num, const OuterBoundary3D &boundary)
{
	vector<Vector3D> points(num);

	for (int i = 0; i < num; i++)
		points[i] = RandomPoint(boundary);

	return points;
}

void RunVoronoi(Tessellation3D *tes, const std::string name)
{
	cout << "Running " << name << "..." << endl;
	tes->Initialise(points, *boundary);
	tessellations[name] = tes;

	cout << "Writing results..." << endl;
	fs::path full_path = fs::path(args.OutputDirectory) / (name + ".out");
	ofstream output;
	output.open(full_path.string());

	for (unsigned i = 0; i < tes->GetPointNo(); i++)
	{
		auto faces = tes->GetCellFaces(i);
		output << "Cell " << namer.GetName(points[i]) << " at " << points[i] << " with " << faces.size() << " faces" << endl;
		// cout << "\tVolume: " << tes.GetVolume(i) << ", Center of Mass: " << namer.GetName(tes.GetCellCM(i)) << tes.GetCellCM(i) << endl;
		for (unsigned j = 0; j < faces.size(); j++)
		{
			auto face = tes->GetFace(faces[j]);

			output << "\tFace F" << faces[j];
			if (face.OtherNeighbor(i).is_initialized())
				output << " neighbor C" << face.OtherNeighbor(i).value().GetCell() + 1;
			output << endl;
			for (auto it = face.vertices.begin(); it != face.vertices.end(); it++)
			{
				output << "\t\t" << namer.GetName(*it) << " " << *it << endl;
			}
			output << endl;
		}
	}
}

static void CompareTessellations()
{
	vector<string> names;

	// Get the keys from tessellations. Remarkably, there's no one function that does this
	pair<string, Tessellation3D *> p;
	BOOST_FOREACH( p, tessellations )
	{
		names.push_back(p.first);
	};

	for (int i = 1; i < tessellations.size(); i++)
		CompareTessellations(names[0], names[i]);
}

#define ERROR(msg) \
	do  \
	{  \
		cout << "****** " << msg << endl; \
		return; \
	} while (false)


static void CompareTessellations(const std::string name1, const std::string name2)
{
	const Tessellation3D *tes1 = tessellations[name1];
	const Tessellation3D *tes2 = tessellations[name2];

	cout << "Comparing " << name1 << " and " << name2 << endl;

	if (tes1->GetPointNo() != tes2->GetPointNo())
		ERROR("Different number of points???");

	for (size_t i = 0; i < tes1->GetPointNo(); i++)
	{
		Vector3D pt1 = tes1->GetMeshPoint(i);
		Vector3D pt2 = tes2->GetMeshPoint(i);
		if (pt1 != pt2)
			ERROR("Different points????");

		CompareCells(tes1, tes2, i);
	}
}

static void CompareCells(const Tessellation3D *tes1, const Tessellation3D *tes2, size_t cellNum)
{
	vector<size_t> faces1 = tes1->GetCellFaces(cellNum);
	vector<size_t> faces2 = tes2->GetCellFaces(cellNum);
	/*if (faces1.size() != faces2.size())
	{
		cout << "***    Wrong number of faces in cell " << cellNum << " " << faces1.size() << " vs. " << faces2.size() << endl;
		return;
	} */


	unordered_set<size_t> faceSet1(faces1.begin(), faces1.end());
	unordered_set<size_t> faceSet2(faces2.begin(), faces2.end());

	// First, look for identical faces
	for (size_t j = 0; j < faces1.size(); j++)
	{
		const Face &face1 = tes1->GetFace(faces1[j]);

		// Find the face in tes2
		bool found = false;
		size_t k = 0;
		for (k = 0; k < faces2.size(); k++)
		{
			const Face &face2 = tes2->GetFace(faces2[k]);
			if (IsSameFace(face1, face2))
			{
				found = true;
				break;
			}
		}

		if (found)
		{
			faceSet1.erase(faces1[j]);
			faceSet2.erase(faces2[k]);
		}
	}

	if (!faceSet1.empty())
	{
		cerr << "***     Can't find matches in C" << cellNum+1 << " for first faces ";
		for (unordered_set<size_t>::iterator it = faceSet1.begin(); it != faceSet1.end(); it++)
			cerr << "F" << *it << " ";
		cerr << endl;
	}
	if (!faceSet2.empty())
	{
		cerr << "***     Can't find matches in C" << cellNum + 1 << " for second faces ";
		for (unordered_set<size_t>::iterator it = faceSet2.begin(); it != faceSet2.end(); it++)
			cerr << "F" << *it << " ";
		cerr << endl;
	}
}


static const double FACE_EPSILON = 1e-5;

// Checks if two faces are the same - this is done by looking at the coordinates - each vertex in
// face1 should appear in face2. We allow close-enough vertices
static bool IsSameFace(const Face &face1, const Face &face2)
{
	if (face1.vertices.size() != face2.vertices.size())
		return false;

	for (vector<VectorRef>::const_iterator it1 = face1.vertices.begin(); it1 != face1.vertices.end(); it1++)
	{
		bool found = false;
		for (vector<VectorRef>::const_iterator it2 = face2.vertices.begin(); it2 != face2.vertices.end(); it2++)
		{
			if (*it1 == *it2)
			{
				found = true;
				break;
			}
			Vector3D diff = **it1 - **it2;
			if (abs(diff) < FACE_EPSILON)
			{
				found = true;
				break;
			}
		}
		if (!found)
			return false;
	}

	return true;
}

static double FaceSimilarity(const Face &face1, const Face &face2)
{
	return -1.0;
}
