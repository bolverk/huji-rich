#include "Lmotion3D.hpp"
#include "../../misc/utils.hpp"

namespace
{
	double estimate_wave_speeds(ComputationalCell3D const& left, ComputationalCell3D const& right,
		EquationOfState const &eos, TracerStickerNames const& tsn)
	{
		const double dl = left.density;
		const double pl = left.pressure;
		const double vl = left.velocity.x;
		const double cl = eos.dp2c(dl, pl, left.tracers, tsn.tracer_names);
		const double dr = right.density;
		const double pr = right.pressure;
		const double vr = right.velocity.x;
		const double cr = eos.dp2c(dr, pr, right.tracers, tsn.tracer_names);
		const double sl = std::min(vl - cl, vr - cr);
		const double sr = std::max(vl + cl, vr + cr);
		const double ss = (pr - pl + dl*vl*(sl - vl) - dr*vr*(sr - vr)) /
			(dl*(sl - vl) - dr*(sr - vr));
		return ss;
	}
}

LMotion3D::LMotion3D(LinearGauss3D const & interp, EquationOfState const & eos,double round_speed,
	double max_v_correction) 
	:interp_(interp), eos_(eos),round_speed_(round_speed),max_v_correction_(max_v_correction){}

void LMotion3D::operator()(const Tessellation3D & tess, const vector<ComputationalCell3D>& cells, double /*time*/,
	TracerStickerNames const & /*tracerstickernames*/, vector<Vector3D>& res) const
{
	size_t N = tess.GetPointNo();
	res.resize(N);
	for (size_t i = 0; i < N; ++i)
		res[i] = cells[i].velocity;
}

void LMotion3D::ApplyFix(Tessellation3D const & tess, vector<ComputationalCell3D> const & cells, double time, double dt,
	vector<Vector3D>& velocities, TracerStickerNames const & tracerstickernames) const
{
	size_t N = tess.GetPointNo();
	velocities.resize(N);
	std::vector<double> TotalArea(N, 0);
	vector<std::pair<ComputationalCell3D, ComputationalCell3D> > edge_values;
	edge_values.resize(tess.GetTotalFacesNumber(), std::pair<ComputationalCell3D, ComputationalCell3D>(cells[0], cells[0]));
	interp_(tess, cells, time, edge_values, tracerstickernames);
	size_t Nfaces = edge_values.size();
	vector<double> ws(Nfaces, 0), face_area(Nfaces, 0);
	std::vector<Vector3D> normals(Nfaces);
	for (size_t j = 0; j < Nfaces; ++j)
	{
		face_area[j] = tess.GetArea(j);
		if (tess.GetFaceNeighbors(j).first < N)
			TotalArea[tess.GetFaceNeighbors(j).first] += face_area[j];
		if (tess.GetFaceNeighbors(j).second < N)
			TotalArea[tess.GetFaceNeighbors(j).second] += face_area[j];
		normals[j] = normalize(tess.GetMeshPoint(tess.GetFaceNeighbors(j).second) -
			tess.GetMeshPoint(tess.GetFaceNeighbors(j).first));
		double par_left = fastabs(edge_values[j].first.velocity - ScalarProd(edge_values[j].first.velocity, normals[j])
			*normals[j]);
		double par_right = fastabs(edge_values[j].second.velocity - ScalarProd(edge_values[j].second.velocity, normals[j])
			*normals[j]);
		edge_values[j].first.velocity.x = ScalarProd(edge_values[j].first.velocity, normals[j]);
		edge_values[j].first.velocity.y = par_left;
		edge_values[j].first.velocity.z = 0;
		edge_values[j].second.velocity.x = ScalarProd(edge_values[j].second.velocity, normals[j]);
		edge_values[j].second.velocity.y = par_right;
		edge_values[j].second.velocity.z = 0;
		ws[j] = estimate_wave_speeds(edge_values[j].first, edge_values[j].second, eos_, tracerstickernames);
	}
	size_t indexX = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaX")) - tracerstickernames.tracer_names.begin());
	size_t indexY = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaY")) - tracerstickernames.tracer_names.begin());
	size_t indexZ = static_cast<size_t>(binary_find(tracerstickernames.tracer_names.begin(), tracerstickernames.tracer_names.end(),
		string("AreaZ")) - tracerstickernames.tracer_names.begin());
	size_t Ntracers = tracerstickernames.tracer_names.size();
	for (size_t i = 0; i < N; ++i)
	{
		face_vec const& faces = tess.GetCellFaces(i);
		double V = tess.GetVolume(i);
		Vector3D CM = V*tess.GetCellCM(i);
		for (size_t j = 0; j < faces.size(); ++j)
		{
			double Vtemp = face_area[faces[j]] * dt*ws[faces[j]];
			if (tess.GetFaceNeighbors(faces[j]).second == i)
				Vtemp *= -1;
			Vector3D CMtemp = tess.FaceCM(faces[j]) + ws[faces[j]] * 0.5*normals[faces[j]] * dt;
			CM += Vtemp*CMtemp;
			V += Vtemp;
		}
		double ratio = fastabs(tess.GetCellCM(i) - tess.GetMeshPoint(i)) / (tess.GetWidth(i));
		double factor = std::min(0.6*ratio*ratio, round_speed_);
		velocities[i] = (CM / V -  tess.GetCellCM(i)) / dt;
		double mag_v = fastabs(velocities[i]);
		Vector3D round_speed = factor * (CM / V - tess.GetMeshPoint(i)) / dt;
		double mag_round = fastabs(round_speed);
		double max_mag_v = 2 * std::max(mag_v, eos_.dp2c(cells[i].density, cells[i].pressure,
			cells[i].tracers, tracerstickernames.tracer_names));
		if (mag_round > max_mag_v)
			round_speed *= max_mag_v / mag_round;
		velocities[i] += round_speed;
		if (indexX < Ntracers && indexY < Ntracers && indexZ < Ntracers)
		{
			double m = 5.5*cells[i].density*tess.GetVolume(i) / (dt*TotalArea[i]);
			Vector3D toadd(m*cells[i].tracers[indexX], m*cells[i].tracers[indexY], m*cells[i].tracers[indexX]);
			double v = fastabs(velocities[i]);
			double t = fastabs(toadd);
			if (t > max_v_correction_*v)
				toadd = toadd*(max_v_correction_*v / t);
			velocities[i] += toadd;
		}
	}
}
