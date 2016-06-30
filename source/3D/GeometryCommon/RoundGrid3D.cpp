#include "RoundGrid3D.hpp"

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
	if (tproc == 0)
		tess->Build(points);
	else
		tess->Build(points, *tproc);
#else
	tess->Build(points);
#endif
	double pi = 3.141592653;
	double eta_ = 0.02, chi_ = 1;
	int N = tess->GetPointNo();
	vector<Vector3D> res(points);
	// Copy the points
	for (size_t i = 0; i < N; ++i)
		res[i] = tess->GetMeshPoint(i);

	for (size_t j = 0; j < NumberIt; ++j)
	{
#ifdef RICH_MPI
		N = tess->GetPointNo();
		res = tess->GetMeshPoints();
		res.resize(static_cast<size_t>(N));
#endif
		for (size_t i = 0; i < N; ++i)
		{
			double R = sqrt(3*tess->GetVolume(i) / (4*pi));
			Vector3D s = tess->GetCellCM(i);
			Vector3D r = tess->GetMeshPoint(i);
			double d = abs(s - r);
			Vector3D dw;
			if (d / eta_ / R < 0.95)
				dw = 0 * s;
			else
				dw = chi_*0.5*(s - r);
			res[i] = tess->GetMeshPoint(i) + dw;
		}
#ifdef RICH_MPI
		if (tproc == 0)
			tess->Build(res);
		else
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
