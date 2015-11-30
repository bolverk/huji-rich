#include "CentroidMotion.hpp"
#include "../../../tessellation/ConvexHull.hpp"
#include "../../../misc/simple_io.hpp"

namespace
{
	vector<size_t> ConvexIndeces(vector<Vector2D> const& points,Vector2D const& meshpoint)
	{
		// Start building the convexhull
		size_t n = points.size();
		vector<double> angles(n);
		for (size_t i = 0; i<n; ++i)
			angles.at(i) = atan2(points.at(i).y - meshpoint.y, points.at(i).x - meshpoint.x);
		return sort_index(angles);
	}

	vector<Vector2D> GetConvexPoints(vector<Vector2D> const& chull, Vector2D const& meshpoint, double R)
	{
		Edge e1, e2;
		vector<Vector2D> res(chull.size());
		for (size_t i = 0; i < res.size(); ++i)
		{
			Vector2D normal = chull[static_cast<size_t>(i)] - meshpoint;
			normal = normal / abs(normal);
			e1.vertices.first = 0.5*(chull[static_cast<size_t>(i)] + meshpoint);
			e1.vertices.first += 50 * R*Vector2D(normal.y, -normal.x);
			e1.vertices.second = e1.vertices.first - 100 * R*Vector2D(normal.y, -normal.x);
			normal = chull[static_cast<size_t>((i+1)%res.size())] - meshpoint;
			normal = normal / abs(normal);
			e2.vertices.first = 0.5*(chull[static_cast<size_t>((i + 1) % res.size())] + meshpoint);
			e2.vertices.first += 50 * R*Vector2D(normal.y, -normal.x);
			e2.vertices.second = e2.vertices.first - 100 * R*Vector2D(normal.y, -normal.x);
			if (!SegmentIntersection(e1, e2, res[i]))
				throw UniversalError("No intersection in centroid motion");
		}
		return res;
	}

	Vector2D GetCM(vector<Vector2D> const& chull,Vector2D const& meshpoint)
	{
		double area = 0;
		Vector2D res;
		for (size_t i = 0; i < chull.size(); ++i)
		{
			const double area_temp = 0.5*std::abs(ScalarProd(chull[i]-meshpoint, zcross(chull[(i+1)%chull.size()]-meshpoint)));
			area += area_temp;
			res += (area_temp / 3.)*(meshpoint+chull[i]+ chull[(i + 1) % chull.size()]);
		}
		return res / area;
	}

	vector<Vector2D> GetOrgChullPoints(Tessellation const& tess,vector<size_t> const& indeces)
	{
		vector<Vector2D> res(indeces.size());
		for (size_t i = 0; i < indeces.size(); ++i)
			res[i] = tess.GetMeshPoint(static_cast<int>(indeces[i]));
		return res;
	}

	vector<Vector2D> GetChullVelocity(Tessellation const & tess, vector<Vector2D> const& orgvel,
		vector<size_t> const& indeces,int point)
	{
		vector<Vector2D> res(indeces.size());
		const size_t N = static_cast<size_t>(tess.GetPointNo());
		for (size_t i = 0; i < indeces.size(); ++i)
		{
			if(indeces[i]<N)
				res[i] = orgvel[indeces[i]];
			else
			{
				if (tess.GetOriginalIndex(static_cast<int>(indeces[i]) == static_cast<int>(indeces[i])))
				{
					Vector2D normal = tess.GetMeshPoint(point) - tess.GetMeshPoint(static_cast<int>(indeces[i]));
					normal = normal / abs(normal);
					Vector2D vel = orgvel[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(indeces[i])))];
					res[i] = vel - 2*normal*ScalarProd(vel,normal);
				}
				else
					res[i] = orgvel[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(indeces[i])))];
			}
		}
		return res;
	}
	
	Vector2D GetCorrectedVelociy(Vector2D const& cm,Vector2D const& meshpoint,Vector2D const& w,
		double dt,double reduce_factor,double R)
	{
		Vector2D temp = cm - meshpoint - dt*w;
		if(abs(temp)<1e-6*R)
			return (cm - meshpoint) / dt;
		const Vector2D newcm = meshpoint + w*dt + std::min(1.0, abs(cm - meshpoint)*reduce_factor/abs(temp))*temp;
		const Vector2D diff = newcm - meshpoint;
		return diff / dt;
	}

	vector<Vector2D> GetCorrectedVelocities(Tessellation const& tess,vector<Vector2D> const& w,double dt,
		double reduce_factor,size_t Niter)
	{
		size_t N = static_cast<size_t>(tess.GetPointNo());
		vector<Vector2D> res(N);
		vector<Vector2D> cur_w(w);
		vector<vector<size_t> > neighbor_indeces(N);
		for (size_t i = 0; i < N; ++i)
		{
			vector<int> neigh = tess.GetNeighbors(static_cast<int>(i));
			for (size_t j = 0; j < neigh.size(); ++j)
				neighbor_indeces[i].push_back(static_cast<size_t>(neigh[j]));
		}
		for (size_t j = 0; j < Niter; ++j)
		{
			for (size_t i = 0; i < N; ++i)
			{
				vector<Vector2D> chull = GetOrgChullPoints(tess, neighbor_indeces[i]);
				vector<size_t> convex_index = ConvexIndeces(chull, tess.GetMeshPoint(static_cast<int>(i)));
				chull = VectorValues(chull,convex_index);
				chull = chull + dt*GetChullVelocity(tess, cur_w,VectorValues(neighbor_indeces[i],convex_index)
					,static_cast<int>(i));
				chull = GetConvexPoints(chull, tess.GetMeshPoint(static_cast<int>(i)) + cur_w[i] * dt,
					tess.GetWidth(static_cast<int>(i)));
				Vector2D cm = GetCM(chull, tess.GetMeshPoint(static_cast<int>(i)) + cur_w[i] * dt);
				res[i] = GetCorrectedVelociy(cm, tess.GetMeshPoint(static_cast<int>(i)), w[i], dt, reduce_factor,
					tess.GetWidth(static_cast<int>(i)));
			}
			cur_w = res;
		}
		// check output
		for (size_t i = 0; i < res.size(); ++i)
		{
			if (tess.GetWidth(static_cast<int>(i))< abs(res[i]) * dt)
				throw UniversalError("Velocity too large in Centroid Motion");
		}
		return res;
	}
}

CentroidMotion::CentroidMotion(double reduction_factor, LinearGaussImproved const& interp, size_t niter) :
	reduce_factor_(reduction_factor),interp_(interp),niter_(niter){}

vector<Vector2D> CentroidMotion::operator()(const Tessellation & tess, const vector<ComputationalCell>& cells, double /*time*/) const
{
	vector<Vector2D> res(static_cast<size_t>(tess.GetPointNo()));
	for (size_t i = 0; i < res.size(); ++i)
	{
		Vector2D density_grad(interp_.GetSlopesUnlimited()[i].first.density, interp_.GetSlopesUnlimited()[i].second.density);
		density_grad = 0.1*density_grad*abs(cells[i].velocity) / abs(density_grad);
		res[i] = cells[i].velocity+density_grad;
	}
	return res;
}

vector<Vector2D> CentroidMotion::ApplyFix(Tessellation const & tess, vector<ComputationalCell> const & /*cells*/, double /*time*/,
	double dt, vector<Vector2D> const & velocities) const
{
	return GetCorrectedVelocities(tess, velocities, dt, reduce_factor_, niter_);
}
