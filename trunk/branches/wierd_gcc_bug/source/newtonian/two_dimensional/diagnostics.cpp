#include "../../misc/simple_io.hpp"
#include "diagnostics.hpp"

namespace
{
	bool VectorSort(Vector2D const& v1,Vector2D const& v2)
	{
		return ((v1.x<v2.x)||((v1.x==v2.x)&&(v1.y<v2.y)));
	}
}

void ConvexHull(vector<Vector2D> &result,Tessellation const* tess,int index)
{
	vector<int> edge_index=tess->GetCellEdges(index);
	result.clear();
	vector<Vector2D> points;
	double R=tess->GetWidth(index);
	points.push_back(tess->GetEdge(edge_index[0]).GetVertex(0));
	points.push_back(tess->GetEdge(edge_index[0]).GetVertex(1));
	// Remove identical points
	for(size_t i=1;i<edge_index.size();++i)
	{
		size_t n=points.size();
		bool samepoint=false;
		for(size_t j=0;j<n;++j)
		{
			if(tess->GetEdge(edge_index[i]).GetVertex(0).distance(points[j])<1e-6*R)
				samepoint=true;
		}
		if(!samepoint)
			points.push_back(tess->GetEdge(edge_index[i]).GetVertex(0));
		samepoint=false;
		for(size_t j=0;j<n;++j)
		{
			if(tess->GetEdge(edge_index[i]).GetVertex(1).distance(points[j])<1e-6*R)
				samepoint=true;
		}
		if(!samepoint)
			points.push_back(tess->GetEdge(edge_index[i]).GetVertex(1));
	}
	// Find the leftmost point
	sort(points.begin(),points.end(),VectorSort);
	// Start building the convexhull
	size_t k=0;
	result.resize(2*points.size());
	double tol=-R*R*1e-9;
	for(size_t i=0;i<points.size();++i)
	{
		while((k>=2)&&(CrossProduct(result[k-1]-result[k-2],points[i]-
			result[k-2])<=tol))
		{
			--k;
		}
		result[k++]=points[i];
	}
	size_t t=k+1;
	for(int i=(int)points.size()-2;i>=0;--i)
	{
		while((k>=t)&&(CrossProduct(result[k-1]-result[k-2],points[i]-
			result[k-2])<=tol))
		{
			--k;
		}
		result[k++]=points[i];
	}
	result.resize(points.size());
}


void DisplayError(UniversalError const& eo,int cycle_number)
{
	cout << eo.GetErrorMessage() << endl;
	cout << "Iteration number = " << cycle_number<< endl;
	for(size_t i=0;i<eo.GetFields().size();++i)
		cout << eo.GetFields()[i] << " = "<< eo.GetValues()[i] << endl;
	throw;
}

namespace {
	double cell_property(hdsim const& sim,int cell_index,string const& property)
	{
		if("generating point x"==property)
			return sim.GetMeshPoint(cell_index).x;
		else if("generating point y"==property)
			return sim.GetMeshPoint(cell_index).y;
		else if("density"==property)
			return sim.GetCell(cell_index).Density;
		else if("pressure"==property)
			return sim.GetCell(cell_index).Pressure;
		else if("velocity x"==property)
			return sim.GetCell(cell_index).Velocity.x;
		else if("velocity y"==property)
			return sim.GetCell(cell_index).Velocity.y;
		else if("volume"==property)
			return sim.GetCellVolume(cell_index);
		else
		  throw UniversalError("Unknown cell property");
	}
}

vector<double> cells_property(hdsim const& sim,
	string const& property)
{
	vector<double> res(sim.GetCellNo(),0);
	for(int i=0;i<sim.GetCellNo();++i)
		res[i] = cell_property(sim,i,property);
	return res;
}

void write_cells_property(hdsim const& sim,
	string const& property,
	string const& fname)
{
	const vector<double> temp = cells_property(sim,property);
	write_vector(temp,fname);
}

string replace_all(string const& base,
	char char_old, char char_new)
{
	string res = base;
	for(int i=0;i<(int)res.size();++i){
		if(res[i]==char_old)
			res[i] = char_new;
	}
	return res;
}

void write_edges_and_neighbors(hdsim const& sim,
	string const& fname)
{
	ofstream f(fname.c_str());
	for(int i=0;i<sim.GetEdgeNo();++i){
		Edge edge = sim.GetEdge(i);
		for(int j=0;j<2;++j)
			f << edge.GetVertex(j).x << " "
			<< edge.GetVertex(j).y << " "
			<< edge.GetNeighbor(j) << " ";
		f << endl;
	}
	f.close();
}

void write_generating_points(hdsim const& sim,
			     string const& fname)
{
  ofstream f(fname.c_str());
  for(int i=0;i<sim.GetCellNo();++i){
    f << cell_property(sim,i,"generating point x") << " "
      << cell_property(sim,i,"generating point y") << endl;
  }
  f.close();
}

void write_array_2d(vector<vector<double> > const& data,
	string const& fname)
{
	ofstream f(fname.c_str());
	for(int i=0;i<(int)data.size();++i){
		for(int j=0;j<(int)data[i].size();++j)
			f << data[i][j] << " ";
		f << endl;
	}
	f.close();
}

Conserved total_conserved(hdsim const& sim)
{
	Conserved res;
	for(int i=0;i<sim.GetCellNo();++i)
		res += Primitive2Conserved
		(sim.GetCell(i),
		sim.GetCellVolume(i));
	return res;
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

