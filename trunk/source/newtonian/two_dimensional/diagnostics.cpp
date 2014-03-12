#include "../../misc/simple_io.hpp"
#include "diagnostics.hpp"

namespace
{
	bool VectorSort(Vector2D const& v1,Vector2D const& v2)
	{
		return ((v1.y<v2.y)||((v1.y==v2.y)&&(v1.x<v2.x)));
	}

	bool AngleSort(Vector2D const& v1,Vector2D const& v2)
	{
		const double tol=1e-8;
		const double a1=atan2(v1.y,v1.x);
		const double a2=atan2(v2.y,v2.x);
		if(a1<a2+tol)
			return true;
		if(a2<a1+tol)
			return false;
		if(abs(v1)<abs(v2))
			return true;
		else
			return false;
	}
}

void ConvexHull(vector<Vector2D> &result,Tessellation const* tess,int index)
{
	vector<int> edge_index=tess->GetCellEdges(index);
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
	// Find the bottom point
	sort(points.begin(),points.end(),VectorSort);
	// Start building the convexhull
	int n=(int)points.size();
	vector<int> indeces(n-1);
	vector<double> angles(n-1);
	for(int i=0;i<n-1;++i)
		angles[i]=atan2(points[i+1].y-points[0].y,points[i+1].x-points[0].x);
	sort_index(angles,indeces);
	result.resize(points.size());
	result[0]=points[0];
	// Check for colinear points
	const double tol=1e-8;
	vector<Vector2D> pfirst,plast;
	for(int i=1;i<n-1;++i)
	{
		if(abs(angles[indeces[i]]-angles[indeces[0]])<tol)
			pfirst.push_back(points[indeces[i]+1]);
		if(abs(angles[indeces[n-i-2]]-angles[indeces[n-2]])<tol)
			plast.push_back(points[indeces[n-i-2]+1]);
	}
	result[0]=points[0];
	if(!pfirst.empty())
	{
		pfirst.insert(pfirst.begin(),points[indeces[0]+1]);
		int N=(int)pfirst.size();
		vector<double> dist(N);
		for(int i=0;i<N;++i)
			dist[i]=abs(pfirst[i]-points[0]);
		vector<int> indeces2(N);
		sort_index(dist,indeces2);
		ReArrangeVector(pfirst,indeces2);
		for(int i=0;i<N;++i)
			result[i+1]=pfirst[i];
	}
	if(!plast.empty())
	{
		plast.insert(plast.begin(),points[indeces[n-2]+1]);
		int N=(int)plast.size();
		vector<double> dist;
		for(int i=0;i<N;++i)
			dist.push_back(abs(plast[i]-points[0]));
		vector<int> indeces2(N);
		sort_index(dist,indeces2);
		ReArrangeVector(plast,indeces2);
		for(int i=0;i<N;++i)
			result[n-1-i]=plast[i];
	}
	int loc1=(int)pfirst.size();
	int loc2=(int)plast.size();
	for(int i=loc1+1;i<n-loc2;++i)
		result[i]=points[indeces[i-1]+1];
}


void DisplayError(UniversalError const& eo)
{
	cout << eo.GetErrorMessage() << endl;
	for(size_t i=0;i<eo.GetFields().size();++i)
		cout << eo.GetFields()[i] << " = "<< eo.GetValues()[i] << endl;
	throw;
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
