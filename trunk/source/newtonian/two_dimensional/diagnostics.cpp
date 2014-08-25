#include "../../misc/simple_io.hpp"
#include "diagnostics.hpp"

void DisplayError(UniversalError const& eo)
{
  cout << eo.GetErrorMessage() << endl;
  for(size_t i=0;i<eo.GetFields().size();++i)
    cout << eo.GetFields()[i] << " = "<< eo.GetValues()[i] << endl;
  throw;
}

void write_error(const string& fname,
		 const UniversalError& eo)
{
  ofstream f(fname.c_str());
  f << eo.GetErrorMessage() << endl;
  for(size_t i=0;eo.GetFields().size();++i)
    f << eo.GetFields()[i] << " = " << eo.GetValues()[i] << endl;
  f.close();
}

Conserved total_conserved(hdsim const& sim)
{
  Conserved res;
  for(int i=0;i<sim.GetCellNo();++i)
    res += Primitive2Conserved
      (sim.GetCell(i),
       sim.GetCellVolume(i));

  #ifdef RICH_MPI
  double total_mass = 0;
  MPI_Allreduce(&res.Mass,&total_mass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  double total_x_momentum = 0;
  MPI_Allreduce(&res.Momentum.x,&total_x_momentum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  double total_y_momentum = 0;
  MPI_Allreduce(&res.Momentum.y,&total_y_momentum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  double total_energy = 0;
  MPI_Allreduce(&res.Energy,&total_energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  res.Mass = total_mass;
  res.Momentum.x = total_x_momentum;
  res.Momentum.y = total_y_momentum;
  res.Energy = total_energy;
  #endif

  return res;
}

vector<Vector2D> ReadVector2DFromFile(string filename)
{
  fstream myFile (filename.c_str(),ios::in | ios::binary);
  if(!myFile.good())
    throw UniversalError("Error opening Vector2D file!!");
  int N;
  myFile.read((char*)&N,sizeof (int));
  vector<Vector2D> res(N);
  for(int i=0;i<N;++i)
    {
      double x,y;
      myFile.read((char*)&x,sizeof(double));
      myFile.read((char*)&y,sizeof(double));
      res[i]=Vector2D(x,y);
    }
  myFile.close();
  return res;
}

void WriteVector2DToFile(vector<Vector2D> const& vec,string filename)
{
  if(vec.empty())
    throw UniversalError("Attempted to write a vector of Vector2D to file with zero length");
  fstream myFile (filename.c_str(),ios::out | ios::binary);
  int n=(int)vec.size();
  myFile.write ((char*)&n,sizeof(int));
  for(int i=0;i<n;++i)
    {
      myFile.write ((char*)&vec[i].x,sizeof(double));
      myFile.write ((char*)&vec[i].y,sizeof(double));
    }
  myFile.close();
}

void BinOutput(string location,hdsim const& sim,Tessellation const& V,bool
	       floatprecision)
{
  fstream myFile (location.c_str(),ios::out | ios::binary);
  int temp=V.GetPointNo();
  myFile.write ((char*)&temp, sizeof (int));
  Vector2D p;
  for(int i=0;i<temp;++i)
    {
      p=V.GetMeshPoint(i);
      if(floatprecision)
	{
	  float x=static_cast<float>(p.x);
	  float y=static_cast<float>(p.y);
	  myFile.write ((char*)&x,sizeof(float));
	  myFile.write ((char*)&y,sizeof(float));
	}
      else
	{
	  myFile.write ((char*)&p.x,sizeof(double));
	  myFile.write ((char*)&p.y,sizeof(double));
	}
    }
  Primitive P;
  for(int i=0;i<temp;++i)
    {
      P=sim.GetCell(i);
      if(floatprecision)
	{
	  float Pressure=static_cast<float>(P.Pressure);
	  myFile.write((char*)&Pressure,sizeof(float));
	}
      else
	myFile.write((char*)&P.Pressure,sizeof(double));
    }
  for(int i=0;i<temp;++i)
    {
      P=sim.GetCell(i);
      if(floatprecision)
	{
	  float Density=static_cast<float>(P.Density);
	  myFile.write((char*)&Density,sizeof(float));
	}
      else
	myFile.write((char*)&P.Density,sizeof(double));
    }
  for(int i=0;i<temp;++i)
    {
      P=sim.GetCell(i);
      if(floatprecision)
	{
	  float xVelocity=static_cast<float>(P.Velocity.x);
	  myFile.write((char*)&xVelocity,sizeof(float));
	}
      else
	myFile.write((char*)&P.Velocity.x,sizeof(double));
    }
  for(int i=0;i<temp;++i)
    {
      P=sim.GetCell(i);
      if(floatprecision)
	{
	  float yVelocity=static_cast<float>(P.Velocity.y);
	  myFile.write((char*)&yVelocity,sizeof(float));
	}
      else
	myFile.write((char*)&P.Velocity.y,sizeof(double));
    }
  // Do the convex hull for each point
  vector<Vector2D> convhull;
  vector<int> nedges(temp);
  for(int i=0;i<temp;++i)
    {
      ConvexHull(convhull,&V,i);
      nedges[i]=(int)convhull.size();
      myFile.write((char*)&nedges[i],sizeof(int));
      for(size_t j=0;j<convhull.size();++j)
	{
	  if(floatprecision)
	    {
	      float x=static_cast<float>(convhull[j].x);
	      float y=static_cast<float>(convhull[j].y);
	      myFile.write((char*)&x,sizeof(float));
	      myFile.write((char*)&y,sizeof(float));
	    }
	  else
	    {
	      myFile.write((char*)&convhull[j].x,sizeof(double));
	      myFile.write((char*)&convhull[j].y,sizeof(double));
	    }
	}
    }
  if(floatprecision)
    {
      float time=static_cast<float>(sim.GetTime());
      myFile.write ((char*)&time,sizeof(float));
    }
  else
    {
      double time=sim.GetTime();
      myFile.write ((char*)&time,sizeof(double));
    }
  const vector<vector<double> > tracers=sim.getTracers();
  int TracerLength=(int)tracers.size();
  myFile.write ((char*)&TracerLength,sizeof(int));
  if(TracerLength>0)
    temp=(int)tracers[0].size();
  else
    temp=0;
  myFile.write ((char*)&temp,sizeof(int));
  float x;
  for(int i=0;i<TracerLength;++i)
    {
      for(int j=0;j<temp;++j)
	{
	  if(floatprecision)
	    {
	      x=static_cast<float>(tracers[i][j]);
	      myFile.write ((char*)&x,sizeof(float));
	    }
	  else
	    {
	      double xx=tracers[i][j];
	      myFile.write ((char*)&xx,sizeof(double));
	    }
	}
    }
  myFile.close();
}