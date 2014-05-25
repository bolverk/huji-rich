#include "CircularRotation.hpp"

CircularRotation::CircularRotation(PointMotion& naive,double Rinner,double Router,
	double w_inner,double w_outer,int Ninner,int Nouter,double t,
	Vector2D const& center):
naive_(naive),inner_radius_(Rinner),outer_radius_(Router),w_inner_(w_inner),
	w_outer_(w_outer),Ninner_(Ninner),Nouter_(Nouter),t0_(t),center_(center),
	init_angles_(vector<double>()),init_R_(vector<double>()),first_time_(true)
{}

Vector2D CircularRotation::CalcVelocity(int index,Tessellation const& tess,
	vector<Primitive> const& cells,double time)
{
	return naive_.CalcVelocity(index,tess,cells,time);
}

vector<Vector2D> CircularRotation::calcAllVelocities(Tessellation const& tess,
	vector<Primitive> const& cells,double time,vector<CustomEvolution*> & cevolve)
{
	if(first_time_)
	{
		first_time_=false;
		init_angles_.resize(Ninner_+Nouter_);
		init_R_.resize(Ninner_+Nouter_);
		for(int i=0;i<Ninner_+Nouter_;++i)
		{
			const Vector2D point=tess.GetMeshPoint(i)-center_;
			init_angles_[i]=atan2(point.y,point.x);
			init_R_[i]=abs(point);
			if(point.y<0)
				init_angles_[i]+=2*M_PI;
		}
	}
	vector<Vector2D> result = naive_.calcAllVelocities(tess,cells,time,cevolve);
	for(int i=0;i<Ninner_+Nouter_;++i)
	{
		const Vector2D point=tess.GetMeshPoint(i)-center_;
		const double r=abs(point);
		const Vector2D phi_hat(-point.y/r,point.x/r);
		const Vector2D r_hat(point.x/r,point.y/r);
		double omega;
		if(abs(r-inner_radius_)<abs(r-outer_radius_))
			omega=w_inner_;
		else
			omega=w_outer_;
		const double vmag=omega*r;
		double p_angle=atan2(point.y,point.x);
		if(point.y<0)
			p_angle+=2*M_PI;
		double time_angle=(time-t0_)*omega+init_angles_[i];
		time_angle-=2*M_PI*floor(time_angle/(2*M_PI));
		double d_angle;
		if(abs(time_angle-p_angle)<0.1)
			d_angle=time_angle-p_angle;
		else
			if(time_angle>p_angle)
				d_angle=(time_angle-2*M_PI)-p_angle;
			else
				d_angle=time_angle-(p_angle-2*M_PI);
		const double R=tess.GetWidth(i);
		if(init_R_[i]>r)
			result[i]=vmag*phi_hat+vmag*sqrt((init_R_[i]-r)/R)*r_hat+
			r*vmag*d_angle*phi_hat/R;
		else
			result[i]=vmag*phi_hat-vmag*sqrt((-init_R_[i]+r)/R)*r_hat+
			r*vmag*d_angle*phi_hat/R;
	}
	return result;
}
