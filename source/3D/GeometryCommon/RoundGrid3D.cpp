#include "RoundGrid3D.hpp"

#include "../../misc/simple_io.hpp"
#include "../../misc/int2str.hpp"
#include <fstream>

vector<Vector3D> RoundGrid3D(vector<Vector3D> const& points, Vector3D const& ll, Vector3D const& ur,
	size_t NumberIt,
#ifdef RICH_MPI
	Tessellation3D const* tproc,
#endif
	Tessellation3D *tess)
{
	Voronoi3D default_tess(ll, ur);
	if (tess == 0)
		tess = &default_tess;
#ifdef RICH_MPI
	tess->Build(points, *tproc);
#else
	tess->Build(points);
#endif
	double eta_ = 0.02, chi_ = 1;
	size_t N = tess->GetPointNo();
	vector<Vector3D> res(points);
	for (size_t j = 0; j < NumberIt; ++j)
	{
#ifdef RICH_MPI
		N = tess->GetPointNo();
		res = tess->GetMeshPoints();
		res.resize(static_cast<size_t>(N));
#endif
		for (size_t i = 0; i < N; ++i)
		{
			double R = tess->GetWidth(i);
			Vector3D s = tess->GetCellCM(i);
			Vector3D r = tess->GetMeshPoint(i);
			double d = fastabs(s - r);
			Vector3D dw;
			if (d / eta_ / R < 0.95)
				dw = 0 * s;
			else
				dw = chi_*0.5*(s - r);
			res[i] = tess->GetMeshPoint(i) + dw;
		}
#ifdef RICH_MPI
		tess->Build(res, *tproc);
#else
		tess->Build(res);
#endif
	}
#ifdef RICH_MPI
	N = tess->GetPointNo();
	res = tess->GetMeshPoints();
	res.resize(static_cast<size_t>(N));
#endif
	return res;
}

#ifdef RICH_MPI
vector<Vector3D> RoundGrid3DSingle(vector<Vector3D> const& points, Vector3D const& ll, Vector3D const& ur,
	size_t NumberIt)
{
	Voronoi3D tess(ll, ur);
	double eta_ = 0.02, chi_ = 1;
	tess.Build(points);
	size_t N = tess.GetPointNo(); //build tess first
	vector<Vector3D> res(points);

	for (size_t j = 0; j < NumberIt; ++j)
	{
		for (size_t i = 0; i < N; ++i)
		{
			double R = tess.GetWidth(i);
			Vector3D s = tess.GetCellCM(i);
			Vector3D r = tess.GetMeshPoint(i);
			double d = fastabs(s - r);
			Vector3D dw;
			if (d / eta_ / R < 0.95)
				dw = 0 * s;
			else
				dw = chi_*0.5*(s - r);
			res[i] = tess.GetMeshPoint(i) + dw;
		}
		if(j<(NumberIt-1))
			tess.Build(res); 
	}
	return res;
}
#endif