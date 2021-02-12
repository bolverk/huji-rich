#include "ANNSelfGravity.hpp"
#include "../../ANN/ANN.h"
#include "../../ANN/kd_tree.h"
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif
#include "../../misc/simple_io.hpp"
#include "../../misc/int2str.hpp"

#ifdef RICH_MPI
namespace
{
	void CallTreeGetSend(ANNkd_tree* tree, Tessellation3D const& tproc, face_vec const& faces, vector<ANNkd_ptr>& nodes,
		double opening)
	{
		point_vec_v facepoints;
		vector<ANNpointArray> annfaces(faces.size());
		vector<size_t> Nfaces(faces.size());
		vector<ANNpoint> normals(faces.size());
		for (size_t j = 0; j < faces.size(); ++j)
		{
			facepoints = VectorValues(tproc.GetFacePoints(), tproc.GetPointsInFace(faces[j]));
			size_t NFace = facepoints.size();
			Nfaces[j] = NFace;
			annfaces[j] = annAllocPts(static_cast<int>(NFace), 3);
#ifdef __INTEL_COMPILER
#pragma omp simd early_exit
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
		tree->GetToSend(annfaces, Nfaces, nodes, opening, normals);
	}
}
#endif
ANNSelfGravity::ANNSelfGravity(double opening, Tessellation3D const* tproc, std::string debug_name) : 
	opening_(opening), tproc_(tproc), d_name_(debug_name) {}

void ANNSelfGravity::operator()(const Tessellation3D & tess, const vector<ComputationalCell3D>& cells,
	const vector<Conserved3D>& /*fluxes*/, const double /*time*/, TracerStickerNames const & /*tracerstickernames*/,
	vector<Vector3D>& acc) const
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
	acc.resize(Norg);
	size_t Nmass = std::max(Norg, static_cast<size_t>(1));
	vector<double> masses(Nmass,0.0);
	std::vector<std::array<double, 6> >  Q(Nmass, std::array<double, 6>{});

	ANNpointArray dpoints = annAllocPts(static_cast<int>(Nmass), 3);
	vector<Vector3D> const& AllCM = tess.GetAllCM();
	vector<double> const& volumes = tess.GetAllVolumes(); 
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
	for (size_t i = 0; i < Norg; ++i)
	{
		dpoints[i][0] = AllCM[i].x;
		dpoints[i][1] = AllCM[i].y;
		dpoints[i][2] = AllCM[i].z;
		masses[i] = volumes[i]*cells[i].density;
	}
	for (size_t i = 0; i < Norg; ++i)
	{
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
		for (size_t j = 0; j < 6; j++)
			Q[i][j] = 0;
	}
#ifdef RICH_MPI
	if (Nmass > Norg)
	{
		dpoints[0][0] = tproc_->GetMeshPoint(rank).x;
		dpoints[0][1] = tproc_->GetMeshPoint(rank).y;
		dpoints[0][2] = tproc_->GetMeshPoint(rank).z;
	}

#endif // RICH_MPI
	ANNkd_tree *atree = new ANNkd_tree(dpoints, masses, Q, static_cast<int>(Nmass), 1, ANN_KD_SL_MIDPT);

#ifdef RICH_MPI
#ifdef timing
	double t1 = MPI_Wtime();
	if (rank == 0)
		std::cout << "Self build time " << t1 - t0 << std::endl;
#endif

	// Get essential tree from other cpu
	size_t Nproc = tproc_->GetPointNo();
	vector<ANNkd_ptr> nodes;
	vector<int> m_size(Nproc, 0);
	vector<double> m_send;

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
			double Qsum = 0;
#ifdef __INTEL_COMPILER
#pragma omp simd reduction(+:Qsum)
#endif
			for (size_t k = 0; k < 6; ++k)
				Qsum += std::fabs(nodes[j]->Q[k]);
			if (Qsum > DBL_MIN)
			{
				for (size_t k = 0; k < 6; ++k)
					m_send.push_back(nodes[j]->Q[k]);
				m_size[i] += 10;
			}
			else
			{
				m_send.push_back(DBL_MAX);
				m_size[i] += 5;
			}
		}
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
	//int max_loc = static_cast<int>(std::max_element(m_size.begin(), m_size.end()) - m_size.begin());
	vector<int> maxes(Nproc);
	MPI_Gather(&Nmaxlocal, 1, MPI_INT, &maxes[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	int Nmaxglobal = 0;
	MPI_Reduce(&Nmaxlocal, &Nmaxglobal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	//	int max_logc = static_cast<int>(std::max_element(maxes.begin(), maxes.end()) - maxes.begin());
	if (rank == 0)
		std::cout << "Maximum send length in gravity is " << Nmaxglobal << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	vector<int> s_disp(Nproc, 0), r_disp(Nproc, 0);
	for (size_t i = 1; i < Nproc; ++i)
		s_disp[i] = s_disp[i - 1] + m_size[i - 1];
	for (size_t i = 1; i < Nproc; ++i)
		r_disp[i] = r_disp[i - 1] + m_rec_size[i - 1];
	vector<double> m_recv(r_disp.back() + m_rec_size.back(), 0);
#ifdef RICH_DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	write_vector(m_send, d_name_ + "m_send_" + int2str(rank) + ".txt");
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
		write_number(static_cast<int>(m_recv.size()), d_name_ + "m_recv_" + int2str(rank) + ".txt");
		throw;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	assert(m_recv.size() % 5 == 0);
	size_t m_recv_size = m_recv.size();
	size_t toadd = 0;
	for (size_t i = 4; i < m_recv_size; i+=5)
	{
		++toadd;
		if (m_recv[i] > (0.1 * DBL_MAX))
			continue;
		else
		{
			i += 5;
			continue;
		}
	}
	dpoints = annAllocPts(static_cast<int>(Nmass + toadd), 3);
	masses.resize(Nmass + toadd);
	Q.resize(Nmass + toadd);
	size_t counter = 0;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t i = 0; i < toadd; ++i)
	{
		masses[Nmass + i] = m_recv[counter];
		dpoints[Nmass + i][0] = m_recv[counter+1];
		dpoints[Nmass + i][1] = m_recv[counter +2];
		dpoints[Nmass + i][2] = m_recv[counter + 3];
		if (m_recv[counter + 4] > (0.1 * DBL_MAX))
		{
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
			for (int j = 0; j < 6; j++)
				Q[Nmass + i][j] = 0;
			counter += 5;
		}
		else
		{
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
			for (int j = 0; j < 6; j++)
				Q[Nmass + i][j] = m_recv[counter + 4 + j];
			counter += 10;
		}
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
	atree = new ANNkd_tree(dpoints, masses, Q, static_cast<int>(Nmass + toadd), 1, ANN_KD_SL_MIDPT);
#ifdef timing
	t0 = MPI_Wtime();
	if (rank == 0)
		std::cout << "Big tree time " << t0 - t1 << std::endl;
#endif
#endif

	// Get acc
	size_t Nbatch = 16;
	ANNpointArray anqpoints = annAllocPts(static_cast<int>(Nbatch), 3);
	ANNpointArray accress = annAllocPts(static_cast<int>(Nbatch), 3);
	std::vector<ANNpoint > qpoints(Nbatch), accpoints(Nbatch);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
	for (size_t i = 0; i < Nbatch; ++i)
	{
		qpoints[i] = anqpoints[i];
		accpoints[i] = accress[i];
	}
	size_t counter2 = 0;
	while (counter2 < Norg)
	{
		size_t Ninner = std::min(Norg - counter2, Nbatch);
		qpoints.resize(Ninner);
		accpoints.resize(Ninner);
#ifdef __INTEL_COMPILER
#pragma omp simd
#endif
		for (size_t j = 0; j < Ninner; ++j)
		{
			qpoints[j][0] = AllCM[counter2+j].x;
			qpoints[j][1] = AllCM[counter2+j].y;
			qpoints[j][2] = AllCM[counter2+j].z;
			accpoints[j][0] = 0;
			accpoints[j][1] = 0;
			accpoints[j][2] = 0;
		}
		counter2+=Ninner;
		atree->GetAcc(qpoints, accpoints, opening_);
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
		for (size_t j = 0; j < Ninner; ++j)
		{
			acc[counter2 - Ninner + j].x = accpoints[j][0];
			acc[counter2 - Ninner + j].y = accpoints[j][1];
			acc[counter2 - Ninner + j].z = accpoints[j][2];
		}
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
