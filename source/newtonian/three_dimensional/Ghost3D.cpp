#include "Ghost3D.hpp"

Ghost3D::~Ghost3D(void) {}

vector<std::pair<size_t, size_t> > Ghost3D::GetOuterFacesIndeces(Tessellation3D const& tess)const
{
	size_t N = tess.GetTotalFacesNumber();
	size_t Npoints = tess.GetPointNo();
	vector<std::pair<size_t, size_t> > res;
	for (size_t i = 0; i < N; ++i)
	{
		if (tess.BoundaryFace(i))
		{
			if(tess.GetFaceNeighbors(i).first>=Npoints)
				res.push_back(std::pair<size_t, size_t>(i, 1));
			else
				res.push_back(std::pair<size_t, size_t>(i, 2));
		}
	}
	return res;
}


boost::container::flat_map<size_t, ComputationalCell3D> RigidWallGenerator3D::operator()(const Tessellation3D& tess,
	const vector<ComputationalCell3D>& cells, double /*time*/, TracerStickerNames const&
	/*tracerstickernames*/) const
{
	boost::container::flat_map<size_t, ComputationalCell3D> res;
	vector<std::pair<size_t, size_t> > ghosts = GetOuterFacesIndeces(tess);
	size_t N = ghosts.size();
	vector<std::pair<size_t, ComputationalCell3D> > temp(N);
	for (size_t i = 0; i < N; ++i)
	{
		size_t ghost_index = 0, real_index = 0;
		if (ghosts[i].second == 1)
		{
			ghost_index = tess.GetFaceNeighbors(ghosts[i].first).first;
			real_index = tess.GetFaceNeighbors(ghosts[i].first).second;
		}
		else
		{
			ghost_index = tess.GetFaceNeighbors(ghosts[i].first).second;
			real_index = tess.GetFaceNeighbors(ghosts[i].first).first;
		}
		Vector3D normal = normalize(tess.GetMeshPoint(ghost_index) - tess.GetMeshPoint(real_index));
		temp[i] = std::pair<size_t, ComputationalCell3D>(ghost_index, cells[real_index]);
		temp[i].second.velocity -= 2*normal*ScalarProd(normal, temp[i].second.velocity);
	}
	res.insert(temp.begin(), temp.end());
	return res;
}

Slope3D RigidWallGenerator3D::GetGhostGradient(const Tessellation3D& /*tess*/, const vector<ComputationalCell3D>& cells,
	const vector<Slope3D>& /*gradients*/, size_t /*ghost_index*/, double /*time*/, size_t /*face_index*/,
	TracerStickerNames const& /*tracerstickernames*/) const
{
	Slope3D res;
	res.xderivative.tracers.resize(cells[0].tracers.size(), 0);
	res.yderivative.tracers.resize(cells[0].tracers.size(), 0);
	res.zderivative.tracers.resize(cells[0].tracers.size(), 0);
	return res;
}

