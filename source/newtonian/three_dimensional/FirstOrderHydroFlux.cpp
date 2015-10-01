#include "FirstOrderHydroFlux.hpp"

namespace
{

  bool approx_equal(double a, double b, double thres=1e-9)
  {
    return thres>std::abs(a-b);
  }

	Vector3D normalize(const Vector3D& v)
	{
		return v/abs(v);
	}

	Vector3D GetParallel(Vector3D const& vel,Vector3D const& normal)
	{
		Vector3D parallel = vel-ScalarProd(vel,normal)*normal;
		const double p=abs(parallel);
		const double v=abs(vel);
		if(approx_equal(v,0))
			return parallel;
		if(p<1e-8*v)
			parallel=parallel/v;
		else
			parallel=parallel/p;
		return parallel;
	}
	Primitive calc_primitive(const ComputationalCell& cc,
		const EquationOfState& eos,
		const Vector3D& normal)
	{
		return Primitive(cc.density,cc.pressure,Vector2D(ScalarProd(cc.velocity,
			normal),abs(cc.velocity-ScalarProd(cc.velocity,normal)*normal)),
			eos.dp2e(cc.density,cc.pressure),eos.dp2c(cc.density,cc.pressure));
	}

	Conserved3D CalcRigidWall(const Tessellation3D& tess,
		const vector<ComputationalCell>& cells,const EquationOfState& eos,
		size_t index,RiemannSolver const& rs)
	{
		Face const& f = tess.GetFace(index);
		const size_t real_index = tess.IsGhostPoint(f.neighbors.first) ? 1 : 0;
		const size_t real_cell = (real_index==1) ? f.neighbors.second : f.neighbors.first;
		const Vector3D normal=normalize(tess.Normal(index));
		Primitive preal=calc_primitive(cells[real_cell],eos, normal);
		const Vector3D parallel = GetParallel(cells[real_cell].velocity,normal);
		ComputationalCell other=cells[real_cell];
		other.velocity=Reflect(cells[real_cell].velocity,normal);
		Primitive pghost = calc_primitive(other,eos,normal);
		const double fv=0;
		Conserved sol;
		if(real_index==0)
			sol = rs(preal,pghost,fv);
		else
			sol = rs(pghost,preal,fv);
		Conserved3D res(sol.Mass,sol.Momentum.x*normal+sol.Momentum.y*parallel,
			sol.Energy);
		if(!other.tracers.empty())
			res.tracers.resize(other.tracers.size(),0);
		return res;
	}

	Conserved3D CalcSingleFluxInBulk(const Tessellation3D& tess,
		const vector<ComputationalCell>& cells,const EquationOfState& eos,
		size_t index,RiemannSolver const& rs,Vector3D const& face_vel)
	{
		Face const& face = tess.GetFace(index);
		const Vector3D normal=normalize(tess.Normal(index));
		const double fv=ScalarProd(face_vel,normal);
		const Conserved sol = rs(calc_primitive(cells[face.neighbors.first],
			eos,normal),calc_primitive(cells[face.neighbors.second],
			eos,normal),fv);
		const ComputationalCell& donor = cells[sol.Mass>0 ? face.neighbors.first :
			face.neighbors.second];
		const Vector3D parallel = GetParallel(donor.velocity,normal);
		Conserved3D res(sol.Mass,sol.Momentum.x*normal+sol.Momentum.y*parallel,
			sol.Energy);
		if(!donor.tracers.empty())
			res.tracers = (sol.Mass>0 ? 1 : -1)*res.mass*donor.tracers;
		return res;
	}

	Conserved3D CalcSingleFlux(const Tessellation3D& tess,
		const vector<ComputationalCell>& cells,const EquationOfState& eos,
		size_t index,RiemannSolver const& rs,Vector3D const& face_vel)
	{
		if(tess.BoundaryFace(index))
			return CalcRigidWall(tess,cells,eos,index,rs);
		return CalcSingleFluxInBulk(tess,cells,eos,index,rs,face_vel);
	}
}

FirstOrderHydroFlux::FirstOrderHydroFlux(RiemannSolver const& rs):rs_(rs){}

FirstOrderHydroFlux::~FirstOrderHydroFlux(void){}

vector<Conserved3D> FirstOrderHydroFlux::operator()(const Tessellation3D& tess,
	const vector<ComputationalCell>& cells,const EquationOfState& eos,
	const vector<Vector3D>& point_velocities) const
{
	vector<Conserved3D> res(tess.GetTotalFacesNumber());
	for(size_t i=0;i<res.size();++i)
	{
		Vector3D fv;
		if(!tess.BoundaryFace(i))
		{
			Face const& f=tess.GetFace(i);
			fv=tess.CalcFaceVelocity(f.neighbors.first,f.neighbors.second,
				point_velocities[f.neighbors.first],point_velocities[f.neighbors.second]);
		}
		res[i]=CalcSingleFlux(tess,cells,eos,i,rs_,fv);
	}
	return res;
}

