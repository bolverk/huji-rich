#include "Reset.hpp"

void ResetOutput(string location,hdsim const& sim)
{
	ResetDump dump;
	fstream myFile (location.c_str(),ios::out | ios::binary);
	int temp=sim.GetCellNo();
	myFile.write ((char*)&temp, sizeof (int));
	Vector2D p;
	for(int i=0;i<temp;i++)
	{
		p=sim.GetMeshPoint(i);
		double x=(p.x);
		double y=(p.y);
		myFile.write ((char*)&x,sizeof(double));
		myFile.write ((char*)&y,sizeof(double));
	}
	Primitive P;
	for(int i=0;i<temp;i++)
	{
		P=sim.GetCell(i);
		double Pressure=(P.Pressure);
		myFile.write((char*)&Pressure,sizeof(double));
		double Density=(P.Density);
		myFile.write((char*)&Density,sizeof(double));
		double xVelocity=(P.Velocity.x);
		myFile.write((char*)&xVelocity,sizeof(double));
		double yVelocity=(P.Velocity.y);
		myFile.write((char*)&yVelocity,sizeof(double));
	}
	double dtemp=sim.GetTime();
	myFile.write ((char*)&dtemp,sizeof(double));
	bool cold=sim.GetColdFlowFlag();
	myFile.write((char*)&cold,sizeof(bool));
	dtemp=sim.GetCfl();
	myFile.write((char*)&dtemp,sizeof(double));
	double a,b;
	sim.GetColdFlowParm(a,b);
	myFile.write((char*)&a,sizeof(double));
	myFile.write((char*)&b,sizeof(double));
	int itemp=sim.GetCycle();
	myFile.write((char*)&itemp,sizeof(int));	
	int n=0;
	vector<vector<double> > tracers=sim.getTracers();
	if(!tracers.empty())
		n=(int)tracers[0].size();
	myFile.write ((char*)&n,sizeof(int));
	if(n==0)
		return;
	if((int)tracers.size()!=temp)
		throw UniversalError("Error in ResetDump, length of mesh points and tracers do not match");
	double x;
	for(int i=0;i<temp;++i)
	{
		for(int j=0;j<n;++j)
		{
			x=tracers[i][j];
			myFile.write ((char*)&x,sizeof(double));
		}
	}
	cold=sim.GetDensityFloorFlag();
	myFile.write((char*)&cold,sizeof(bool));
	double dtemp2;
	sim.GetDensityFloorParm(dtemp,dtemp2);
	myFile.write((char*)&dtemp,sizeof(double));
	myFile.write((char*)&dtemp2,sizeof(double));
	myFile.close();
}

void ResetRead(string location,ResetDump &dump,EquationOfState const* eos)
{
	fstream myFile (location.c_str(),ios::in | ios::binary);
	if(!myFile.good())
		throw UniversalError("Error opening reset file!!");
	int N;
	myFile.read((char*)&N,sizeof (int));
	// Resize the vectors
	dump.snapshot.cells.resize(N);
	dump.snapshot.mesh_points.resize(N);
	dump.tracers.clear();
	// Read the data
	Vector2D cortemp;
	double x,y,d,p;
	Primitive ptemp;
	for(int i=0;i<N;++i)
	{
		myFile.read((char*)&x,sizeof(double));
		myFile.read((char*)&y,sizeof(double));
		cortemp.Set(x,y);
		dump.snapshot.mesh_points[i]=cortemp;
	}
	for(int i=0;i<N;++i)
	{
		myFile.read((char*)&p,sizeof(double));
		myFile.read((char*)&d,sizeof(double));
		myFile.read((char*)&cortemp.x,sizeof(double));
		myFile.read((char*)&cortemp.y,sizeof(double));
		dump.snapshot.cells[i]=CalcPrimitive(d,p,cortemp,*eos);
	}
	myFile.read((char*)&dump.time,sizeof(double));
	myFile.read((char*)&dump.coldflows,sizeof(bool));
	myFile.read((char*)&dump.cfl,sizeof(double));
	myFile.read((char*)&dump.a,sizeof(double));
	myFile.read((char*)&dump.b,sizeof(double));
	myFile.read((char*)&dump.cycle,sizeof(int));
	int n;
	myFile.read((char*)&n,sizeof(int));
	if(n==0)
		return;
	dump.tracers.resize(N);
	for(int i=0;i<N;++i)
	{
		dump.tracers[i].resize(n);
		for(int j=0;j<n;++j)
		{
			myFile.read((char*)&x,sizeof(double));
			dump.tracers[i][j]=x;
		}
	}
	myFile.read((char*)&dump.densityfloor,sizeof(bool));
	myFile.read((char*)&dump.densitymin,sizeof(double));
	myFile.read((char*)&dump.pressuremin,sizeof(double));
	myFile.close();
}
