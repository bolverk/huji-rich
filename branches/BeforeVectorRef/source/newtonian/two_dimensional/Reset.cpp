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
#ifdef RICH_MPI
	int temp2=sim.GetProcTessellation().GetPointNo();
	myFile.write ((char*)&temp2, sizeof (int));
	for(int i=0;i<temp2;i++)
	{
		p=sim.GetProcTessellation().GetMeshPoint(i);
		double x=(p.x);
		double y=(p.y);
		myFile.write ((char*)&x,sizeof(double));
		myFile.write ((char*)&y,sizeof(double));
	}
#endif

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
	char cold= sim.GetColdFlowFlag() ? '1' : '0';
	myFile.write((char*)&cold,sizeof(char));
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
	cold = static_cast<char>(sim.GetDensityFloorFlag());
	myFile.write((char*)&cold, sizeof(char));
	double dtemp2;
	sim.GetDensityFloorParm(dtemp, dtemp2);
	myFile.write((char*)&dtemp, sizeof(double));
	myFile.write((char*)&dtemp2, sizeof(double));
	for (int i = 0; i<temp; ++i)
	{
		unsigned int ctemp = (unsigned int)sim.custom_evolution_indices[static_cast<size_t>(i)];
		myFile.write((char*)&ctemp, sizeof(unsigned int));
	}
	if (n == 0)
	{
		myFile.close();
		return;
	}
	if((int)tracers.size()!=temp)
		throw UniversalError("Error in ResetDump, length of mesh points and tracers do not match");
	double x;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<temp;++j)
		{
		  x=tracers[static_cast<size_t>(j)][static_cast<size_t>(i)];
			myFile.write ((char*)&x,sizeof(double));
		}
	}
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
	dump.snapshot.cells.resize(static_cast<size_t>(N));
	dump.snapshot.mesh_points.resize(static_cast<size_t>(N));
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
		dump.snapshot.mesh_points[static_cast<size_t>(i)]=cortemp;
	}

#ifdef RICH_MPI
	int temp2;
	myFile.read((char*)&temp2,sizeof (int));
	dump.procmesh.resize(temp2);
	for(int i=0;i<temp2;i++)
	{
		myFile.read((char*)&x,sizeof(double));
		myFile.read((char*)&y,sizeof(double));
		cortemp.Set(x,y);
		dump.procmesh[i]=cortemp;
	}
#endif

	for(int i=0;i<N;++i)
	{
		myFile.read((char*)&p,sizeof(double));
		myFile.read((char*)&d,sizeof(double));
		myFile.read((char*)&cortemp.x,sizeof(double));
		myFile.read((char*)&cortemp.y,sizeof(double));
		dump.snapshot.cells[static_cast<size_t>(i)]=CalcPrimitive(d,p,cortemp,*eos);
	}
	myFile.read((char*)&dump.time,sizeof(double));
	char ctemp;
	myFile.read((char*)&ctemp,sizeof(char));
	if (ctemp == '1')
		dump.coldflows = true;
	else
		dump.coldflows = false;
	myFile.read((char*)&dump.cfl,sizeof(double));
	myFile.read((char*)&dump.a,sizeof(double));
	myFile.read((char*)&dump.b,sizeof(double));
	myFile.read((char*)&dump.cycle,sizeof(int));
	int n;
	myFile.read((char*)&n,sizeof(int));
	myFile.read((char*)&ctemp, sizeof(char));
	dump.densityfloor = static_cast<bool>(ctemp);
	myFile.read((char*)&dump.densitymin, sizeof(double));
	myFile.read((char*)&dump.pressuremin, sizeof(double));
	dump.cevolve.resize(static_cast<size_t>(N));
	for (int i = 0; i<N; ++i)
		myFile.read((char*)&dump.cevolve[static_cast<size_t>(i)], sizeof(unsigned int));
	if (n == 0)
	{
		myFile.close();
		return;
	}
	dump.tracers.resize(static_cast<size_t>(N));
	for(int i=0;i<n;++i)
	{
	  dump.tracers[static_cast<size_t>(i)].resize(static_cast<size_t>(n));
		for(int j=0;j<N;++j)
		{
			myFile.read((char*)&x,sizeof(double));
			dump.tracers[static_cast<size_t>(j)][static_cast<size_t>(i)]=x;
		}
	}
	myFile.close();
}
