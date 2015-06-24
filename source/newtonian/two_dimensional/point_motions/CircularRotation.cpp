#include "CircularRotation.hpp"

KeplerianOmega::KeplerianOmega(double Mass, double RigidMin, double RigidMax) : mass_(Mass), RigidMin_(RigidMin),
RigidMax_(RigidMax){}

double KeplerianOmega::CalcOmega(Vector2D const& point) const
{
	double R = std::min(std::max(abs(point), RigidMin_), RigidMax_);
	return sqrt(mass_ / (R*R*R));
}


CircularRotation::CircularRotation(OmegaFunction const& omega) : evencall_(false), omega_(omega){}

Vector2D CircularRotation::CalcVelocity(int index, Tessellation const& tess,
	vector<Primitive> const& /*cells*/, double /*time*/)
{
	double w = omega_.CalcOmega(tess.GetMeshPoint(index));
	Vector2D const& point = tess.GetMeshPoint(index);
	return Vector2D(-w*point.y, w*point.x);
}


void CircularRotation::ApplyFix(Tessellation const& tess, vector<Primitive> const& /*cells*/, double /*time*/,
	vector<CustomEvolution*> &/*cevolve*/, const vector<vector<double> >& /*tracers*/, double dt, vector < Vector2D >
	& velocities)
{
	size_t N = static_cast<size_t> (tess.GetPointNo());
	for (size_t i = 0; i<N; ++i)
	{
		const Vector2D point = tess.GetMeshPoint(static_cast<int>(i));
		double omega = omega_.CalcOmega(tess.GetMeshPoint((int)i));
		if (!evencall_)
			velocities[i] = Vector2D(-point.x + point.x*cos(omega*dt) - point.y*sin(omega*dt), -point.y + point.y*cos(dt*omega) +
			point.x*sin(dt*omega)) / dt;
		else
			velocities[i] = Vector2D(-2 * point.y*sin(omega*dt*0.5), 2 * point.x*sin(omega*dt*0.5)) / dt;
	}
	evencall_ = !evencall_;
}