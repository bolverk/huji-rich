#include "CentroidMotion.hpp"
#include "../../../tessellation/ConvexHull.hpp"
#include "../../../misc/simple_io.hpp"
#ifdef RICH_MPI
#include <mpi.h>
#endif

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
			e1.vertices.first += 500 * R*Vector2D(normal.y, -normal.x);
			e1.vertices.second = e1.vertices.first - 1000 * R*Vector2D(normal.y, -normal.x);
			normal = chull[static_cast<size_t>((i+1)%res.size())] - meshpoint;
			normal = normal / abs(normal);
			e2.vertices.first = 0.5*(chull[static_cast<size_t>((i + 1) % res.size())] + meshpoint);
			e2.vertices.first += 500 * R*Vector2D(normal.y, -normal.x);
			e2.vertices.second = e2.vertices.first - 1000 * R*Vector2D(normal.y, -normal.x);
			if (!SegmentIntersection(e1, e2, res[i]))
			{
				UniversalError eo("No intersection in centroid motion");
#ifdef RICH_MPI
				int rank;
				MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				eo.AddEntry("Rank", rank);
#endif
				eo.AddEntry("mesh.x", meshpoint.x);
				eo.AddEntry("mesh.y", meshpoint.y);
				eo.AddEntry("R", R);
				for (size_t k = 0; k < chull.size(); ++k)
				{
					eo.AddEntry("chull x", chull[k].x);
					eo.AddEntry("chull y", chull[k].y);
				}
				eo.AddEntry("e1.vertices.first.x", e1.vertices.first.x);
				eo.AddEntry("e1.vertices.first.y", e1.vertices.first.y);
				eo.AddEntry("e1.vertices.second.x", e1.vertices.second.x);
				eo.AddEntry("e1.vertices.second.y", e1.vertices.second.y);
				eo.AddEntry("e2.vertices.first.x", e2.vertices.first.x);
				eo.AddEntry("e2.vertices.first.y", e2.vertices.first.y);
				eo.AddEntry("e2.vertices.second.x", e2.vertices.second.x);
				eo.AddEntry("e2.vertices.second.y", e2.vertices.second.y);
				throw eo;
			}
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
		vector<size_t> const& indeces,
#ifdef RICH_MPI
	int /*point*/
#else
		int point
#endif
		)
	{
		vector<Vector2D> res(indeces.size());
		const size_t N = static_cast<size_t>(tess.GetPointNo());
		for (size_t i = 0; i < indeces.size(); ++i)
		{
			if(indeces[i]<N)
				res[i] = orgvel[indeces[i]];
			else
			{
#ifdef RICH_MPI
				res[i] = orgvel[indeces[i]];
#else
				if (tess.GetOriginalIndex(static_cast<int>(indeces[i]) == static_cast<int>(indeces[i])))
				{
					Vector2D normal = tess.GetMeshPoint(point) - tess.GetMeshPoint(static_cast<int>(indeces[i]));
					normal = normal / abs(normal);
					Vector2D vel = orgvel[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(indeces[i])))];
					res[i] = vel - 2*normal*ScalarProd(vel,normal);
				}
				else

					res[i] = orgvel[static_cast<size_t>(tess.GetOriginalIndex(static_cast<int>(indeces[i])))];
#endif
			}
		}
		return res;
	}
	
	Vector2D GetCorrectedVelociy(Vector2D const& cm,Vector2D const& meshpoint,Vector2D const& w,
		double dt,double reduce_factor,double R,EquationOfState const& eos,ComputationalCell const& cell)
	{
		Vector2D temp = cm - meshpoint - dt*w;
		if(abs(temp)<1e-6*R)
			return (cm - meshpoint) / dt;
		const double distance_reduce = std::min(abs(temp) / R,1.0);
		const double toadd_velocity = reduce_factor*distance_reduce * std::max(abs(w), eos.dp2c(cell.density, cell.pressure));
		return (w + temp*(toadd_velocity / abs(temp)));
	}

	bool ShouldCalc(vector<string> const& toignore, ComputationalCell const& cell,TracerStickerNames const&
		tracerstickernames)
	{
		for (size_t i = 0; i < toignore.size(); ++i)
		{
			vector<string>::const_iterator it = binary_find(tracerstickernames.sticker_names.begin(),
				tracerstickernames.sticker_names.end(), toignore[i]);
			if (it != tracerstickernames.sticker_names.end() && cell.stickers[static_cast<size_t>(
				it - tracerstickernames.sticker_names.begin())])
				return false;
		}
		return true;
	}

	vector<Vector2D> GetCorrectedVelocities(Tessellation const& tess,vector<Vector2D> const& w,double dt,
		double reduce_factor,size_t Niter,vector<ComputationalCell> const& cells,EquationOfState const& eos,
		vector<string> const& toignore, TracerStickerNames const&	tracerstickernames)
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
				if (!ShouldCalc(toignore, cells[i],tracerstickernames))
					continue;
				vector<Vector2D> chull = GetOrgChullPoints(tess, neighbor_indeces[i]);
				vector<size_t> convex_index = ConvexIndeces(chull, tess.GetMeshPoint(static_cast<int>(i)));
				chull = VectorValues(chull,convex_index);
				chull = chull + dt*GetChullVelocity(tess, cur_w,VectorValues(neighbor_indeces[i],convex_index)
					,static_cast<int>(i));
				chull = GetConvexPoints(chull, tess.GetMeshPoint(static_cast<int>(i)) + cur_w[i] * dt,
					tess.GetWidth(static_cast<int>(i)));
				Vector2D cm = GetCM(chull, tess.GetMeshPoint(static_cast<int>(i)) + cur_w[i] * dt);
				res[i] = GetCorrectedVelociy(cm, tess.GetMeshPoint(static_cast<int>(i)), w[i], dt, reduce_factor,
					tess.GetWidth(static_cast<int>(i)),eos,cells[i]);
			}
			cur_w = res;
		}
		// check output
		for (size_t i = 0; i < res.size(); ++i)
		{
			if (tess.GetWidth(static_cast<int>(i)) < abs(res[i]) * dt && abs(cells[i].velocity)*2 < abs(res[i]))
				throw UniversalError("Velocity too large in Centroid Motion");
		}
		return res;
	}
}

CentroidMotion::CentroidMotion(PointMotion const& bpm,double reduction_factor, EquationOfState const& eos, size_t niter,
	vector<string> toignore) :
	bpm_(bpm),reduce_factor_(reduction_factor),eos_(eos),niter_(niter),toignore_(toignore){}

vector<Vector2D> CentroidMotion::operator()(const Tessellation & tess, const vector<ComputationalCell>& cells,
	double time, TracerStickerNames const& tracerstickernames) const
{
	vector<Vector2D> res=bpm_(tess, cells, time, tracerstickernames);
	return res;
}

vector<Vector2D> CentroidMotion::ApplyFix(Tessellation const & tess, vector<ComputationalCell> const & cells, double /*time*/,
	double dt, vector<Vector2D> const & velocities, TracerStickerNames const& tracerstickernames) const
{
	return GetCorrectedVelocities(tess, velocities, dt, reduce_factor_, niter_,cells,eos_,toignore_,
		tracerstickernames);
}
