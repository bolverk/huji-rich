#include "OpticalDepthCalc.hpp"
#include "../../ANN/ANN.h"
#include "../..//ANN/kd_tree.h"
#ifdef RICH_MPI
#include "../..//mpi/mpi_commands.hpp"
#endif
#include "../../misc/simple_io.hpp"
#include "../../misc/int2str.hpp"

namespace
{
	double CumSumDepth(std::vector<std::pair<double, double> >& surface_density, Vector3D const& cm,
		double R)
	{
		size_t N = surface_density.size();
		assert(N > 0);
		if (N == 1)
			return 0.5 * surface_density[0].second * R;
		double sum = 0;
		size_t index = 0;
		for (size_t i = 0; i < N; ++i)
		{
			if (std::abs(surface_density[i].first - cm.z) < 1e-7 * R)
			{
				index = i;
				break;
			}
			surface_density[i].first -= cm.z;
			surface_density[i].second += sum;
			sum = surface_density[i].second;
		}
		if (index == N - 1) // when cell is the heighest
		{
			double result = surface_density[N - 1].second * R * 0.5;
			return result;
		}
		if (index == 0) // when cell is the lowest
		{
			double result = surface_density[0].second * R * 0.5;
			return result;
		}
		sum = 0;
		double density_self = 0.5 * surface_density[index].second;
		for (size_t i = N - 1; i > index; i--)
		{
			surface_density[i].first -= cm.z;
			surface_density[i].second += sum;
			sum = surface_density[i].second;
		}
		if (surface_density.at(index - 1).second > surface_density.at(index + 1).second)
		{
			double total_density = surface_density[index - 1].second + density_self;
			for (int i = static_cast<int>(index) - 1; i >= 0; --i)
			{
				if (surface_density[i].second < 0.5 * total_density || i == 0)
				{
					double result = -total_density * (surface_density[i].first - 0.5 * R);
					return result;
				}
			}
		}
		else
		{
			double total_density = surface_density[index + 1].second + density_self;
			for (size_t i = index + 1; i < N; ++i)
			{
				if (surface_density[i].second < 0.5 * total_density || i == N - 1)
				{
					double result = total_density * (surface_density[i].first + 0.5 * R);
					return result;
				}
			}
		}
		assert(false);
	}

#ifdef RICH_MPI
	void CallTreeGetSend(ANNkd_tree* tree, Tessellation3D const& tproc, face_vec const& faces, vector<ANNkd_ptr>& nodes,
		double opening)
	{
		point_vec_v facepoints;
		std::vector<ANNpointArray> annfaces(faces.size());
		std::vector<size_t> Nfaces(faces.size());
		std::vector<ANNpoint> normals(faces.size());
		for (size_t j = 0; j < faces.size(); ++j)
		{
			facepoints = VectorValues(tproc.GetFacePoints(), tproc.GetPointsInFace(faces[j]));
			size_t NFace = facepoints.size();
			Nfaces[j] = NFace;
			annfaces[j] = annAllocPts(static_cast<int>(NFace), 3);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
			for (size_t i = 0; i < NFace; ++i)
			{
				annfaces[j][i][0] = facepoints[i].x;
				annfaces[j][i][1] = facepoints[i].y;
				annfaces[j][i][2] = facepoints[i].z;
			}
			Vector3D vnorm = normalize(tproc.GetMeshPoint(tproc.GetFaceNeighbors(faces[j]).second) -
				tproc.GetMeshPoint(tproc.GetFaceNeighbors(faces[j]).first));
			normals[j][0] = vnorm.x;
			normals[j][1] = vnorm.y;
			normals[j][2] = vnorm.z;
		}
		tree->GetToSendOpticalDepth(annfaces, Nfaces, nodes, opening);
	}
#endif

}
OpticalDepthCalc::OpticalDepthCalc(double opening, Tessellation3D const* tproc, std::string debug_name) : 
	opening_(opening), tproc_(tproc), d_name_(debug_name) {}

void OpticalDepthCalc::operator()(const Tessellation3D& tess, const vector<ComputationalCell3D>& cells, std::vector<double>& res) const
{
#ifdef RICH_MPI
	assert(tproc_ != 0);
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#ifdef timing
	double t0 = MPI_Wtime();
#endif
#endif
	size_t Norg = tess.GetPointNo();
	res.resize(Norg);
	size_t Nmass = std::max(Norg, static_cast<size_t>(1));
	vector<double> masses(Nmass, 0.0);

	ANNpointArray dpoints = annAllocPts(static_cast<int>(Nmass), 3);
	vector<Vector3D> const& AllCM = tess.GetAllCM();
	vector<double> const& volumes = tess.GetAllVolumes();
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t i = 0; i < Norg; ++i)
	{
		dpoints[i][0] = AllCM[i].x;
		dpoints[i][1] = AllCM[i].y;
		dpoints[i][2] = AllCM[i].z;
		masses[i] = volumes[i] * cells[i].density;
	}
#ifdef RICH_MPI
	if (Nmass > Norg)
	{
		dpoints[0][0] = tproc_->GetMeshPoint(rank).x;
		dpoints[0][1] = tproc_->GetMeshPoint(rank).y;
		dpoints[0][2] = tproc_->GetMeshPoint(rank).z;
	}
#endif // RICH_MPI

	ANNkd_tree* atree = new ANNkd_tree(dpoints, masses, static_cast<int>(Nmass), 1);

#ifdef RICH_MPI
#ifdef timing
	double t1 = MPI_Wtime();
	if (rank == 0)
		std::cout << "Self build time " << t1 - t0 << std::endl;
#endif

	// Get essential tree from other cpu
	size_t Nproc = tproc_->GetPointNo();
	std::vector<ANNkd_ptr> nodes;
	std::vector<int> m_size(Nproc, 0);
	std::vector<double> m_send;

	for (size_t i = 0; i < Nproc; ++i)
	{
		if (i == static_cast<size_t>(rank))
			continue;
		nodes.clear();
		face_vec const& faces = tproc_->GetCellFaces(i);
		CallTreeGetSend(atree, *tproc_, faces, nodes, opening_);
		for (size_t j = 0; j < nodes.size(); ++j)
		{
			m_send.push_back(nodes[j]->mass);
			m_send.push_back(nodes[j]->CM[0]);
			m_send.push_back(nodes[j]->CM[1]);
			m_send.push_back(nodes[j]->CM[2]);
		}
		m_size[i] = static_cast<int>(nodes.size() * 4);
		assert(m_size[i] > 0);
	}
	delete atree;
	annClose();
#ifdef timing
	t0 = MPI_Wtime();
	if (rank == 0)
		std::cout << "Self send calc time " << t0 - t1 << std::endl;
#endif

	// send/recv data
	vector<int> m_rec_size(Nproc);
	MPI_Alltoall(&m_size[0], 1, MPI_INT, &m_rec_size[0], 1, MPI_INT, MPI_COMM_WORLD);

	int Nmaxlocal = *std::max_element(m_size.begin(), m_size.end());
	//	int max_loc = static_cast<int>(std::max_element(m_size.begin(), m_size.end()) - m_size.begin());
	vector<int> maxes(Nproc);
	MPI_Gather(&Nmaxlocal, 1, MPI_INT, &maxes[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	int Nmaxglobal = 0;
	MPI_Reduce(&Nmaxlocal, &Nmaxglobal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	//	int max_logc = static_cast<int>(std::max_element(maxes.begin(), maxes.end()) - maxes.begin());
	if (rank == 0)
		std::cout << "Maximum send length in opticaldepth is " << Nmaxglobal << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	vector<int> s_disp(Nproc, 0), r_disp(Nproc, 0);
	for (size_t i = 1; i < Nproc; ++i)
		s_disp[i] = s_disp[i - 1] + m_size[i - 1];
	for (size_t i = 1; i < Nproc; ++i)
		r_disp[i] = r_disp[i - 1] + m_rec_size[i - 1];
	vector<double> m_recv(r_disp.back() + m_rec_size.back(), 0);
#ifdef RICH_DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	write_vector(m_send, d_name_ +  "m_send_" + int2str(rank) + ".txt");
	write_vector(m_size, d_name_ + "m_size_" + int2str(rank) + ".txt");
	write_vector(s_disp, d_name_ + "s_disp_" + int2str(rank) + ".txt");
	write_vector(m_rec_size, d_name_ + "m_rec_size_" + int2str(rank) + ".txt");
	write_vector(r_disp, d_name_ + "r_disp_" + int2str(rank) + ".txt");
	write_number(m_recv.size(), d_name_ + "m_recv_" + int2str(rank) + ".txt");
#endif
	int error = MPI_Alltoallv(&m_send[0], &m_size[0], &s_disp[0], MPI_DOUBLE, &m_recv[0], &m_rec_size[0], &r_disp[0], MPI_DOUBLE, MPI_COMM_WORLD);
	if (error != 0)
	{
		write_vector(m_send, d_name_ + "m_send_" + int2str(rank) + ".txt");
		write_vector(m_size, d_name_ + "m_size_" + int2str(rank) + ".txt");
		write_vector(s_disp, d_name_ + "s_disp_" + int2str(rank) + ".txt");
		write_vector(m_rec_size, d_name_ + "m_rec_size_" + int2str(rank) + ".txt");
		write_vector(r_disp, d_name_ + "r_disp_" + int2str(rank) + ".txt");
		write_number(static_cast<int>(m_recv.size()), d_name_ + "m_recv_" + int2str(static_cast<int>(rank)) + ".txt");
		throw;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	assert(m_recv.size() % 4 == 0);
	size_t toadd = m_recv.size() / 4;
	dpoints = annAllocPts(static_cast<int>(Nmass + toadd), 3);
	masses.resize(Nmass + toadd);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t i = 0; i < toadd; ++i)
	{
		masses[Nmass + i] = m_recv[i * 4];
		dpoints[Nmass + i][0] = m_recv[i * 4 + 1];
		dpoints[Nmass + i][1] = m_recv[i * 4 + 2];
		dpoints[Nmass + i][2] = m_recv[i * 4 + 3];
	}
#ifdef timing
	t1 = MPI_Wtime();
	if (rank == 0)
		std::cout << "Comm time " << t1 - t0 << std::endl;
#endif

	//rebuild tree	
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t i = 0; i < Norg; ++i)
	{
		dpoints[i][0] = AllCM[i].x;
		dpoints[i][1] = AllCM[i].y;
		dpoints[i][2] = AllCM[i].z;
	}
#ifdef RICH_MPI
	if (Nmass > Norg)
	{
		dpoints[0][0] = tproc_->GetMeshPoint(rank).x;
		dpoints[0][1] = tproc_->GetMeshPoint(rank).y;
		dpoints[0][2] = tproc_->GetMeshPoint(rank).z;
	}
#endif // RICH_MPI
	atree = new ANNkd_tree(dpoints, masses, static_cast<int>(Nmass + toadd), 1);
#ifdef timing
	t0 = MPI_Wtime();
	if (rank == 0)
		std::cout << "Big tree time " << t0 - t1 << std::endl;
#endif
#endif

	// Get depth
	ANNpoint anqpoints;
	std::vector<std::pair<double, double> > tau_res; // distance, dtau
	for (size_t i = 0; i < Norg; ++i)
	{
		tau_res.clear();
		anqpoints[0] = AllCM[i].x;
		anqpoints[1] = AllCM[i].y;
		anqpoints[2] = AllCM[i].z;
		atree->GetOpticalDepth(anqpoints, tau_res, opening_);
		std::sort(tau_res.begin(), tau_res.end());
		res[i] = CumSumDepth(tau_res, AllCM[i], tess.GetWidth(i));
	}
#ifdef RICH_MPI
#ifdef timing
	t1 = MPI_Wtime();
	if (rank == 0)
		std::cout << "Walk time " << t1 - t0 << std::endl;
#endif
#endif
	// Cleanup
	delete atree;
	annClose();
}

