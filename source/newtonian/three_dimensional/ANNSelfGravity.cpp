#include "ANNSelfGravity.hpp"
#include "../../ANN/ANN.h"
#include "../../ANN/kd_tree.h"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif

#ifdef RICH_MPI
namespace
{
	void CallTreeGetSend(ANNkd_tree* tree, Tessellation3D const& tproc, std::vector<size_t>const& faces, vector<ANNkd_ptr>& nodes,
		double opening)
	{
		vector<Vector3D> facepoints;
		vector<ANNpointArray> annfaces(faces.size());
		vector<size_t> Nfaces(faces.size());
		vector<ANNpoint> normals(faces.size());
		for (size_t j = 0; j < faces.size(); ++j)
		{
			facepoints = VectorValues(tproc.GetFacePoints(), tproc.GetPointsInFace(faces[j]));
			size_t NFace = facepoints.size();
			Nfaces[j] = NFace;
			annfaces[j] = annAllocPts(static_cast<int>(NFace), 3);
			for (size_t i = 0; i < NFace; ++i)
			{
				annfaces[j][i][0] = facepoints[i].x;
				annfaces[j][i][1] = facepoints[i].y;
				annfaces[j][i][2] = facepoints[i].z;
			}
			Vector3D vnorm = normalize(tproc.GetMeshPoint(tproc.GetFaceNeighbors(faces[j]).second) -
				tproc.GetMeshPoint(tproc.GetFaceNeighbors(faces[j]).first));
			normals[j] = annAllocPt(3);
			normals[j][0] = vnorm.x;
			normals[j][1] = vnorm.y;
			normals[j][2] = vnorm.z;
		}
		tree->GetToSend(annfaces, Nfaces, nodes, opening, normals);
		for (size_t j = 0; j < faces.size(); ++j)
		{
			annDeallocPts(annfaces[j]);
			annDeallocPt(normals[j]);
		}
		
	}
}
#endif
ANNSelfGravity::ANNSelfGravity(double opening,Tessellation3D const* tproc) : opening_(opening),tproc_(tproc){}

void ANNSelfGravity::operator()(const Tessellation3D & tess, const vector<ComputationalCell3D>& cells,
	const vector<Conserved3D>& /*fluxes*/, const double /*time*/, TracerStickerNames const & /*tracerstickernames*/,
	vector<Vector3D>& acc) const
{
	size_t Norg = tess.GetPointNo();
	acc.resize(Norg);
	vector<double> masses(Norg);
	std::vector<boost::array<double, 6> >  Q(Norg);

	ANNpointArray dpoints = annAllocPts(static_cast<int>(Norg), 3);
	for (size_t i = 0; i < Norg; ++i)
	{
		Vector3D const& vec = tess.GetCellCM(i);
		dpoints[i][0] = vec.x;
		dpoints[i][1] = vec.y;
		dpoints[i][2] = vec.z;
		masses[i] = tess.GetVolume(i)*cells[i].density;
		Q[i].assign(0);
	}
	ANNkd_tree *atree = new ANNkd_tree(dpoints, masses,Q, static_cast<int>(Norg), 3, 1, ANN_KD_SL_MIDPT);

#ifdef RICH_MPI
	// Get essential tree from other cpu
	assert(tproc_ != 0);
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	size_t Nproc = tproc_->GetPointNo();
	vector<vector<double> > send_masses,send_Q;
	vector<vector<Vector3D> > send_CM;
	vector<int> totalkwith;
	vector<ANNkd_ptr> nodes;
	vector<double> mass_temp;
	vector<Vector3D> cm_temp;
	vector<double> Q_temp;
	/*vector<int> m_size(Nproc, 0);
	vector<double> m_send,CM_send,Q_send;
	size_t counter = 0;*/
	for (size_t i = 0; i < Nproc; ++i)
	{
		if (i == static_cast<size_t>(rank))
			continue;
		nodes.clear();
		vector<size_t> const& faces = tproc_->GetCellFaces(i);
		CallTreeGetSend(atree, *tproc_, faces, nodes, opening_);
		cm_temp.clear();
		mass_temp.clear();
		Q_temp.clear();
		for (size_t j = 0; j < nodes.size(); ++j)
		{
			cm_temp.push_back(Vector3D(nodes[j]->CM[0], nodes[j]->CM[1], nodes[j]->CM[2]));
			mass_temp.push_back(nodes[j]->mass);
			for (size_t k = 0; k < 6; ++k)
				Q_temp.push_back(nodes[j]->Q[k]);
			/*m_send.push_back(nodes[j]->mass);
			CM_send.push_back(nodes[j]->CM[0]);
			CM_send.push_back(nodes[j]->CM[1]);
			CM_send.push_back(nodes[j]->CM[2]);*/
		}
		send_CM.push_back(cm_temp);
		send_masses.push_back(mass_temp);
		send_Q.push_back(Q_temp);
		//m_size[i] = static_cast<int>(send_masses.size());
		totalkwith.push_back(static_cast<int>(i));
	}
	// send/recv data
	//vector<int> m_rec(Norg);
	//MPI_Alltoall(&m_size[0],1,MPI_INT, &m_rec[0], 1, MPI_INT,MPI_COMM_WORLD);
	send_masses=MPI_exchange_data(totalkwith, send_masses);
	send_Q = MPI_exchange_data(totalkwith, send_Q);
	send_CM = MPI_exchange_data(totalkwith, send_CM,Vector3D());
	//rebuild tree
	annDeallocPts(dpoints);
	delete atree;
	annClose();
	size_t toadd = 0;
	for (size_t i = 0; i < send_masses.size(); ++i)
		toadd += send_masses[i].size();
	dpoints = annAllocPts(static_cast<int>(Norg+toadd), 3);
	masses.resize(Norg + toadd);
	Q.resize(Norg + toadd);
	for (size_t i = 0; i < Norg; ++i)
	{
		Vector3D const& vec = tess.GetCellCM(i);
		dpoints[i][0] = vec.x;
		dpoints[i][1] = vec.y;
		dpoints[i][2] = vec.z;
	}
	size_t counter = Norg;
	for (size_t i = 0; i < send_masses.size(); ++i)
	{
		for (size_t j = 0; j < send_masses[i].size(); ++j)
		{
			dpoints[counter][0] = send_CM[i][j].x;
			dpoints[counter][1] = send_CM[i][j].y;
			dpoints[counter][2] = send_CM[i][j].z;
			masses[counter] = send_masses[i][j];
			for (size_t k = 0; k < 6; ++k)
				Q[counter][k] = send_Q[i][j * 6 + k];
			++counter;
		}
	}
	atree = new ANNkd_tree(dpoints, masses,Q, static_cast<int>(Norg+toadd), 3, 1, ANN_KD_SL_MIDPT);
#endif

	// Get acc

	size_t Nbatch = 8;
	ANNpointArray anqpoints = annAllocPts(static_cast<int>(Nbatch), 3);
	ANNpointArray accress = annAllocPts(static_cast<int>(Nbatch), 3);
	std::vector<ANNpoint> qpoints(Nbatch), accpoints(Nbatch);
	for (size_t i = 0; i < Nbatch; ++i)
	{
		qpoints[i] = anqpoints[i];
		accpoints[i] = accress[i];
	}
	size_t counter = 0;
	while(counter<Norg)
	{
		size_t Ninner = std::min(Norg - counter, Nbatch);
		qpoints.resize(Ninner);
		accpoints.resize(Ninner);
		for (int j = 0; j < Ninner; ++j)
		{
			Vector3D const& vec = tess.GetCellCM(counter);
			qpoints[j][0] = vec.x;
			qpoints[j][1] = vec.y;
			qpoints[j][2] = vec.z;
			accpoints[j][0] = 0;
			accpoints[j][1] = 0;
			accpoints[j][2] = 0;		
			++counter;
		}
		atree->GetAcc(qpoints, accpoints, opening_);
		for (int j = 0; j < Ninner; ++j)
		{
			acc[counter - Ninner + j].x = accpoints[j][0];
			acc[counter - Ninner + j].y = accpoints[j][1];
			acc[counter - Ninner + j].z = accpoints[j][2];
		}
	}
	annDeallocPts(anqpoints);
	annDeallocPts(accress);	
/*
	ANNpoint anqpoint, accres;
	anqpoint = annAllocPt(3);
	accres = annAllocPt(3);
	for (size_t i = 0; i < Norg; ++i)
	{
		Vector3D const& vec = tess.GetCellCM(i);
		anqpoint[0] = vec.x;
		anqpoint[1] = vec.y;
		anqpoint[2] = vec.z;
		accres[0] = 0;
		accres[1] = 0;
		accres[2] = 0;
		atree->GetAcc(anqpoint, accres, opening_);
		acc[i].x = accres[0];
		acc[i].y = accres[1];
		acc[i].z = accres[2];
	}
	annDeallocPt(anqpoint);
	annDeallocPt(accres);
	*/
	// Cleanup
	annDeallocPts(dpoints);
	delete atree;
	annClose();

	// Fix the effect of close neighbors and add smoothing length
	vector<size_t> neigh;
	vector<double> volumes(Norg);
	for (size_t i = 0; i < Norg; ++i)
		volumes[i] = tess.GetVolume(i);
#ifdef RICH_MPI
	MPI_exchange_data2(tess, volumes, true);
#endif
	vector<double> smoothlength(volumes.size());
	for (size_t i = 0; i < volumes.size(); ++i)
		smoothlength[i] = std::pow(volumes[i], 0.33333333);

	for (size_t i = 0; i < Norg; ++i)
	{
		Vector3D myCM = tess.GetCellCM(i);
		tess.GetNeighbors(i, neigh);
		size_t Nneigh = neigh.size();
		for (size_t j = 0; j < Nneigh; ++j)
		{
			if (neigh[j] >= Norg && tess.IsPointOutsideBox(neigh[j]))
				continue;
			double otherM = volumes[neigh[j]] * cells[neigh[j]].density;
			double distance = abs(myCM - tess.GetCellCM(neigh[j]));
			double new_distance = smoothlength[i] * smoothlength[i] + smoothlength[neigh[j]] * smoothlength[neigh[j]];
			acc[i] += otherM*(myCM - tess.GetCellCM(neigh[j])) *(1.0/(distance*distance*distance) - 
				1.0/(new_distance*sqrt(new_distance)));
		}
	}
}
