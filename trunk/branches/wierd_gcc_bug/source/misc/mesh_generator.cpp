#include "mesh_generator.hpp"

vector<Vector2D> RandPointsRmax(int PointNum,double Rmin,double Rmax,
	double xc,double yc)
{
	double ran[2];
	int counter=0;
	vector<Vector2D> res;
	Vector2D point;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	while(counter<PointNum)
	{
		ran[0]=dist(generator);
		ran[1]=dist(generator);
		point.x=((Rmax-Rmin)*ran[0]+Rmin)*cos(ran[1]*2*M_PI)+xc;
		point.y=((Rmax-Rmin)*ran[0]+Rmin)*sin(ran[1]*2*M_PI)+yc;
		res.push_back(point);
		++counter;
	}	
	return res;
}


vector<Vector2D> RandSquare(int PointNum,double xl,double xr,double yd,
	double yu)
{
	double ran[2];
	double xc=xr-xl;
	double yc=yu-yd;
	vector<Vector2D> res;
	Vector2D point;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	for(int i=0;i<PointNum;++i)
	{
		ran[0]=dist(generator);
		ran[1]=dist(generator);
		point.x=ran[0]*xc+xl;
		point.y=ran[1]*yc+yd;
		res.push_back(point);
	}	
	return res;
}


vector<Vector2D> RandPointsR(int PointNum,double xl,double xr,double yd,
	double yu,double minR)
{
	double ran[2];
	int counter=0;
	double maxR= sqrt(2.0)*max((xr-xl)/2,(yu-yd)/2);
	double xc=(xl+xr)*0.5;
	double yc=(yu+yd)*0.5;
	vector<Vector2D> res;
	Vector2D point;
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	while(counter<PointNum)
	{
		ran[0]=dist(generator);
		ran[1]=dist(generator);
		point.x=maxR*ran[0]*cos(ran[1]*2*M_PI)+xc;
		point.y=maxR*ran[0]*sin(ran[1]*2*M_PI)+yc;
		if(abs(point)>minR)
		{
			if(point.x>xl&&point.x<xr)
				if(point.y>yd&&point.y<yu)
				{
					res.push_back(point);
					++counter;
				}
		}
	}	
	return res;
}



vector<Vector2D> SquarePertubed(int nx,int ny,double sidex,double sidey,
	double mag)
{
	base_generator_type generator;
	boost::random::uniform_real_distribution<> dist;
	double ran[2];
	vector<Vector2D> res(nx*ny);
	double widthx = sidex/(double)nx;
	double widthy = sidey/(double)ny;
	Vector2D point;
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			ran[0]=2*dist(generator)-1;
			ran[1]=2*dist(generator)-1;
			point.x = ((double)i+0.5+mag*ran[0])*widthx-sidex/2;
			point.y = ((double)j+0.5+mag*ran[1])*widthy-sidey/2;
			res[i*ny+j] = point;
		}
	}
	return res;
}


vector<Vector2D> CirclePointsRmax_1(int PointNum,double Rmin,double Rmax,
	double xc,double yc,double xmax,double ymax,double xmin,double ymin)
{
	double Nphid=sqrt(PointNum*2*M_PI/log(Rmax/Rmin));
	int Nr=(int)(Nphid*log(Rmax/Rmin)/(2*M_PI));
	int Nphi=(int)floor(Nphid+0.5);
	double dphi=2*M_PI/Nphi;
	Vector2D pos;
	vector<Vector2D> res;
	double r;
	for(int i=0;i<Nr;++i)
	{
		r=Rmin*exp(2*M_PI*i/Nphi);
		for(int j=0;j<Nphi;++j)
		{
			pos.Set(r*cos(dphi*j)+xc,r*sin(dphi*j)+yc);
			if(pos.x<xmax&&pos.x>xmin)
				if(pos.y<ymax&&pos.y>ymin)
					res.push_back(pos);
		}
	}
	return res;
}


vector<Vector2D> CirclePointsRmax_2(int PointNum,double Rmin,double Rmax,
	double xc,double yc,double xmax,double ymax,double xmin,double ymin)
{
	double x=sqrt(Rmin/Rmax);
	int Nr=int(sqrt(2*PointNum*(1-x)/(1+x)/M_PI));
	double A=4*pow(1/sqrt(Rmin)-1/sqrt(Rmax),2)/(Nr*Nr);
	double dphi;
	int Nphi;
	Vector2D pos;
	vector<Vector2D> res;
	double r;
	for(int i=0;i<Nr;++i)
	{
		r=4*Rmin/pow(2-sqrt(A*Rmin)*i,2);
		dphi=sqrt(A*r);
		Nphi=int(2*M_PI/dphi);
		dphi=double(2*M_PI/Nphi);
		for(int j=0;j<Nphi;++j)
		{
			pos.Set(r*cos(dphi*j)+xc,r*sin(dphi*j)+yc);
			if(pos.x<xmax&&pos.x>xmin)
				if(pos.y<ymax&&pos.y>ymin)
					res.push_back(pos);
		}
	}
	return res;
}

vector<Vector2D> CirclePointsRmax(int PointNum,double Rmin,double Rmax,
	double xc,double yc)
{
	double A=sqrt(M_PI*(Rmax*Rmax-Rmin*Rmin)/PointNum);
	int Nr=int((Rmax-Rmin)/A);
	double dr=(Rmax-Rmin)/Nr;
	Vector2D pos;
	vector<Vector2D> res;
	for(int i=0;i<Nr;++i)
	{
		double r=Rmin+i*dr;
		int Nphi=int(2*M_PI*r/A);
		double dphi=A/r;
		for(int j=0;j<Nphi;++j)
		{
			pos.Set(r*cos(dphi*j)+xc,r*sin(dphi*j)+yc);
			res.push_back(pos);
		}
	}
	return res;
}

vector<Vector2D> SquareMesh(int nx,int ny,double sidex,double sidey,bool centered)
{
	vector<Vector2D> res(nx*ny);
	double widthx = sidex/(double)nx;
	double widthy = sidey/(double)ny;
	Vector2D point;
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			point.x = ((double)i+0.5)*widthx-sidex/2;
			point.y = ((double)j+0.5)*widthy-sidey/2;
			if(!centered)
			{
				point.x+=sidex/2;
				point.y+=sidey/2;
			}
			res[i*ny+j] = point;
		}
	}
	return res;
}

vector<Vector2D> cartesian_mesh(int nx, int ny,
				Vector2D const& lower_left,
				Vector2D const& upper_right)
{
  assert(nx>0);
  assert(ny>0);
  assert(upper_right.x>lower_left.x);
  assert(upper_right.y>lower_left.y);

  vector<Vector2D> res(nx*ny);
  const double dx = (upper_right.x-lower_left.x)/(double)nx;
  const double dy = (upper_right.y-lower_left.y)/(double)ny;
  for(int i=0;i<nx;++i)
    {
      for(int j=0;j<ny;++j)
	{
	  Vector2D point(((double)i+0.5)*dx+lower_left.x,
			 ((double)j+0.5)*dy+lower_left.y);
	  res[i*ny+j] = point;
	}
    }
  return res;
}

vector<Vector2D> Circle(int np,double R,double xc,double yc)
{
	double dphi=2*M_PI/np;
	Vector2D res2;
	vector<Vector2D> res;
	for(int i=0;i<np;++i)
	{
		res2.Set(R*cos(dphi*i)+xc,R*sin(dphi*i)+yc);
		res.push_back(res2);
	}
	return res;
}

vector<Vector2D> Line(int PointNum,double xmin,double xmax,double ymin,double ymax)
{
	double dy=ymax-ymin,dx=xmax-xmin;
	double angle=atan2(dy,dx);
	double length=sqrt(dy*dy+dx*dx);
	double dl=length/PointNum;

	vector<Vector2D> res;
	Vector2D temp;

	for(int i=0;i<PointNum;++i)
	{
		temp.Set(xmin+dl*(i+0.5)*cos(angle),ymin+dl*(i+0.5)*sin(angle));
		res.push_back(temp);
	}
	return res;
}

vector<Vector2D> CirclePointsRmax_a(int PointNum,double Rmin,double Rmax,
	double xc,double yc,double xmax,double ymax,double xmin,double ymin,
	double alpha)
{
	double N0=sqrt(PointNum*4*M_PI*(alpha+1)/(pow(Rmax,2*(alpha+1))-
		pow(Rmin,2*(alpha+1))));
	int Nr=int(floor((pow(Rmax,alpha+1)-pow(Rmin,alpha+1))*N0/(2*M_PI*(alpha+1))));
	double dphi;
	int Nphi;
	Vector2D pos;
	vector<Vector2D> res;
	double r;
	for(int i=0;i<Nr;++i)
	{
		r=pow(2*M_PI*i*(alpha+1)/N0+pow(Rmin,alpha+1),1.0/(alpha+1));
		Nphi=int(floor(N0*pow(r,1+alpha)+1.5));
		dphi=2*M_PI/Nphi;
		for(int j=0;j<Nphi;++j)
		{
			pos.Set(r*cos(dphi*j)+xc,r*sin(dphi*j)+yc);
			if(pos.x<xmax&&pos.x>xmin)
				if(pos.y<ymax&&pos.y>ymin)
					res.push_back(pos);
		}
	}
	return res;
}

