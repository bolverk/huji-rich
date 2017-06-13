#include "TreeSelfGravity.hpp"
#include <cfloat>
extern "C"
{
#include "rebound.h"
#include "tools.h"
#include "output.h"
}
#ifdef RICH_MPI
#include "../../mpi/mpi_commands.hpp"
#endif
TreeSelfGravity::TreeSelfGravity(size_t nx, size_t ny, size_t nz) :nx_(nx), ny_(ny), nz_(nz) {}

void TreeSelfGravity::operator()(const Tessellation3D & tess, const vector<ComputationalCell3D>& cells,
	const vector<Conserved3D>& /*fluxes*/, const double /*time*/, TracerStickerNames const & /*tracerstickernames*/,
	vector<Vector3D>& acc) const
{
	size_t Norg = tess.GetPointNo();
	std::pair<Vector3D, Vector3D> box = tess.GetBoxCoordinates();
	Vector3D CM = 0.5*(box.first + box.second);
	Vector3D diff = box.second - box.first;
	double boxsize = std::max(std::max(diff.x, diff.y), diff.z);

	struct reb_simulation* const r = reb_create_simulation();
	// Setup constants
	r->integrator = reb_simulation::REB_INTEGRATOR_NONE;
	r->gravity = reb_simulation::REB_GRAVITY_TREE;
	r->boundary = reb_simulation::REB_BOUNDARY_NONE;
	r->opening_angle2 = 0.5;          // This constant determines the accuracy of the tree code gravity estimate.
	r->G = 1;
	r->softening = 0;         // Gravitational softening length
	r->dt = DBL_MIN * 100;         // Timestep
								   // Setup root boxes for gravity tree.
								   // Here, we use 2x2=4 root boxes (each with length 'boxsize')
								   // This allows you to use up to 4 MPI nodes.
	reb_configure_box(r, boxsize, static_cast<int>(nx_), static_cast<int>(ny_), static_cast<int>(nz_));


#ifdef RICH_MPI
	int ws = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	reb_mpi_init(r);
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	/*vector<uint32_t> hash(Norg);
	for (size_t i = 0; i < Norg; ++i)
		hash[i] = static_cast<uint32_t>(ws*i + rank);
	// create send data
	Vector3D v3temp;
	vector<vector<double> > send_masses(8), send_x(8), send_y(8), send_z(8);
	vector<vector<uint32_t> > send_h(8);
	// create send data
	for (size_t i = 0; i < Norg; ++i)
	{
		v3temp = tess.GetMeshPoint(i);
		v3temp -= CM;
		int ii = ((int)floor((v3temp.x +boxsize) / boxsize) + 2) % 2;
		int j = ((int)floor((v3temp.y + boxsize) / boxsize) + 2) % 2;
		int k = ((int)floor((v3temp.z + boxsize) / boxsize) + 2) % 2;
		size_t index = static_cast<size_t>((k*2 + j)*2 + ii);
		send_masses[index].push_back(tess.GetVolume(i)*cells[i].density);
		send_x[index].push_back(v3temp.x);
		send_y[index].push_back(v3temp.y);
		send_z[index].push_back(v3temp.z);
		send_h[index].push_back(hash[i]);
	}
	vector<MPI_Request> req(40);
	// send data
	for (int i = 0; i < 8; ++i)
	{
		double dtemp;
		uint32_t utemp;
		if (send_x[static_cast<size_t>(i)].empty())
		{
			MPI_Isend(&dtemp, 0, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5)]);
			MPI_Isend(&dtemp, 0, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5 + 1)]);
			MPI_Isend(&dtemp, 0, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5 + 2)]);
			MPI_Isend(&dtemp, 0, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5 + 3)]);
			MPI_Isend(&utemp, 0, MPI_UINT32_T, i, 4, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5 + 4)]);
		}
		else
		{
			MPI_Isend(&send_x[static_cast<size_t>(i)][0], static_cast<int>(send_x[static_cast<size_t>(i)].size()), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5)]);
			MPI_Isend(&send_y[static_cast<size_t>(i)][0], static_cast<int>(send_y[static_cast<size_t>(i)].size()), MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5 + 1)]);
			MPI_Isend(&send_z[static_cast<size_t>(i)][0], static_cast<int>(send_z[static_cast<size_t>(i)].size()), MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5 + 2)]);
			MPI_Isend(&send_masses[static_cast<size_t>(i)][0], static_cast<int>(send_masses[static_cast<size_t>(i)].size()), MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5 + 3)]);
			MPI_Isend(&send_h[static_cast<size_t>(i)][0], static_cast<int>(send_h[static_cast<size_t>(i)].size()), MPI_UINT32_T, i, 4, MPI_COMM_WORLD, &req[static_cast<size_t>(i * 5 + 4)]);
		}
	}
	// recv data and run tree and then send data
	if (rank < 8)
	{
		MPI_Status status;
		int count;
		vector<double> d_temp;
		vector<double> m_recv, x_recv, y_recv, z_recv;
		vector<uint32_t> u_recv,u_temp;
		for (size_t i = 0; i < static_cast<size_t>(ws * 5); ++i)
		{
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if (status.MPI_TAG < 4)
			{
				MPI_Get_count(&status, MPI_DOUBLE, &count);
				d_temp.resize(static_cast<size_t>(count));
				MPI_Recv(&d_temp[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (status.MPI_TAG == 0)
					x_recv.insert(x_recv.end(), d_temp.begin(), d_temp.end());
				if (status.MPI_TAG == 1)
					y_recv.insert(y_recv.end(), d_temp.begin(), d_temp.end());
				if (status.MPI_TAG == 2)
					z_recv.insert(z_recv.end(), d_temp.begin(), d_temp.end());
				if (status.MPI_TAG == 3)
					m_recv.insert(m_recv.end(), d_temp.begin(), d_temp.end());
			}
			else
			{
				assert(status.MPI_TAG == 4);
				MPI_Get_count(&status, MPI_UINT32_T, &count);
				u_temp.resize(static_cast<size_t>(count));
				MPI_Recv(&u_temp[0], count, MPI_UINT32_T, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				u_recv.insert(u_recv.end(), u_temp.begin(), u_temp.end());
			}
		}
		struct reb_particle pt = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
		size_t Nadd = m_recv.size();
		assert(m_recv.size() == x_recv.size());
		assert(m_recv.size() == y_recv.size());
		assert(m_recv.size() == z_recv.size());
		assert(m_recv.size() == u_recv.size());
		for (size_t i = 0; i < Nadd; ++i)
		{
			pt.x = x_recv[i];
			pt.y = y_recv[i];
			pt.z = z_recv[i];
			pt.m = m_recv[i];
			pt.hash = u_recv[i];
			reb_add(r, pt);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		reb_integrate(r, r->dt);
		assert(static_cast<size_t>(r->N) == Nadd);
		send_h.clear();
		send_x.clear();
		send_y.clear();
		send_z.clear();
		send_h.resize(static_cast<size_t>(ws));
		send_x.resize(static_cast<size_t>(ws));
		send_y.resize(static_cast<size_t>(ws));
		send_z.resize(static_cast<size_t>(ws));
		for (size_t i = 0; i < Nadd; ++i)
		{
			size_t proc_index = static_cast<size_t>(r->particles[i].hash) % static_cast<size_t>(ws);
			uint32_t p_index = static_cast<uint32_t>(static_cast<size_t>(r->particles[i].hash)/ static_cast<size_t>(ws));
			send_h[proc_index].push_back(p_index);
			send_x[proc_index].push_back(r->particles[i].ax);
			send_y[proc_index].push_back(r->particles[i].ay);
			send_z[proc_index].push_back(r->particles[i].az);
		}
		req.resize(static_cast<size_t>(ws * 4));
		for (size_t i = 0; i < static_cast<size_t>(ws); ++i)
		{
			double dtemp;
			uint32_t utemp;
			if (send_h[i].empty())
			{
				MPI_Isend(&utemp,0, MPI_UINT32_T, static_cast<int>(i), 0, MPI_COMM_WORLD, &req[i * 4]);
				MPI_Isend(&dtemp, 0, MPI_DOUBLE, static_cast<int>(i), 1, MPI_COMM_WORLD, &req[i * 4 + 1]);
				MPI_Isend(&dtemp, 0, MPI_DOUBLE, static_cast<int>(i), 2, MPI_COMM_WORLD, &req[i * 4 + 2]);
				MPI_Isend(&dtemp, 0, MPI_DOUBLE, static_cast<int>(i), 3, MPI_COMM_WORLD, &req[i * 4 + 3]);
			}
			else
			{
				MPI_Isend(&send_h[i][0], static_cast<int>(send_h[i].size()), MPI_UINT32_T, static_cast<int>(i), 0, MPI_COMM_WORLD, &req[i * 4]);
				MPI_Isend(&send_x[i][0], static_cast<int>(send_x[i].size()), MPI_DOUBLE, static_cast<int>(i), 1, MPI_COMM_WORLD, &req[i * 4 + 1]);
				MPI_Isend(&send_y[i][0], static_cast<int>(send_y[i].size()), MPI_DOUBLE, static_cast<int>(i), 2, MPI_COMM_WORLD, &req[i * 4 + 2]);
				MPI_Isend(&send_z[i][0], static_cast<int>(send_z[i].size()), MPI_DOUBLE, static_cast<int>(i), 3, MPI_COMM_WORLD, &req[i * 4 + 3]);
			}
		}
	}
	// recv the data from the tree code
	vector<vector<double> > recv_x(8), recv_y(8), recv_z(8);
	vector<double> d_temp;
	vector<uint32_t> u_temp;
	double dtemp;
	uint32_t utemp;
	vector<vector<uint32_t> > recv_h(8);
	MPI_Status status;
	int count;
	acc.resize(Norg);
	hash.clear();
	for (size_t i = 0; i < 8*4; ++i)
	{
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if (status.MPI_TAG == 0)
		{
			MPI_Get_count(&status, MPI_UINT32_T, &count);
			if (count > 0)
			{
				u_temp.resize(static_cast<size_t>(count));
				MPI_Recv(&u_temp[0], count, MPI_UINT32_T, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				recv_h[static_cast<size_t>(status.MPI_SOURCE)] = u_temp;
			}
			else
				MPI_Recv(&utemp, 0, MPI_UINT32_T, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Get_count(&status, MPI_DOUBLE, &count);
			if (count > 0)
			{
				d_temp.resize(static_cast<size_t>(count));
				MPI_Recv(&d_temp[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if(status.MPI_TAG==1)
					recv_x[static_cast<size_t>(status.MPI_SOURCE)] = d_temp;
				if (status.MPI_TAG == 2)
					recv_y[static_cast<size_t>(status.MPI_SOURCE)] = d_temp;
				if (status.MPI_TAG == 3)
					recv_z[static_cast<size_t>(status.MPI_SOURCE)] = d_temp;
			}
			else
				MPI_Recv(&dtemp, 0, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	// Convert recv data to acc
	for (size_t i = 0; i < 8; ++i)
	{
		size_t Nvector = recv_h[i].size();
		for (size_t j = 0; j < Nvector; ++j)
		{
			size_t uindex = static_cast<size_t> (recv_h[i][j]);
			acc[uindex].Set(recv_x[i][j], recv_y[i][j], recv_z[i][j]);
		}
	}
	*/
	struct reb_particle pt = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	for (size_t i = 0; i < Norg; ++i)
	{
		Vector3D const& point = tess.GetMeshPoint(i);
		pt.x = point.x - CM.x;
		pt.y = point.y - CM.y;
		pt.z = point.z - CM.z;
		pt.vx = 0;
		pt.vy = 0;
		pt.vz = 0;
		pt.m = cells[i].density*tess.GetVolume(i);
		pt.hash = static_cast<uint32_t>(ws*i + rank);
		reb_add(r, pt);
	}
	// Start the integration
	reb_integrate(r, r->dt);

	// send back the particles
	size_t WS = static_cast<size_t>(ws);
	vector<vector<uint32_t> > indeces(WS);
	vector<vector<double> > ax(WS), ay(WS), az(WS);
	vector<int> talk_num(WS, 0);
	size_t Nrebound = static_cast<size_t>(r->N);
	for (size_t i = 0; i < Nrebound; ++i)
	{
		size_t proc_index = static_cast<size_t>(r->particles[i].hash) % WS;
		uint32_t local_index = static_cast<uint32_t>(static_cast<size_t>(r->particles[i].hash) / WS);
		indeces[proc_index].push_back(local_index);
		ax[proc_index].push_back(r->particles[i].ax);
		ay[proc_index].push_back(r->particles[i].ay);
		az[proc_index].push_back(r->particles[i].az);
	}
	for (size_t i = 0; i < WS; ++i)
		if (!indeces[i].empty())
			talk_num[i] = 1;
	vector<int> scounts(WS, 1);
	int ntalkwithme = 0;
	MPI_Reduce_scatter(&talk_num[0], &ntalkwithme, &scounts[0], MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Request req;
	for (size_t i = 0; i < WS; ++i)
	{
		int itemp=0;
		if (!indeces[i].empty())
			MPI_Isend(&itemp, 1, MPI_INT, static_cast<int>(i), 0, MPI_COMM_WORLD, &req);
	}
	vector<int> talkwithme;
	for (int i = 0; i < ntalkwithme; ++i)
	{
		MPI_Status status;
		int itemp=0;
		MPI_Recv(&itemp, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		talkwithme.push_back(status.MPI_SOURCE);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	sort(talkwithme.begin(), talkwithme.end());
	for (size_t i = 0; i < WS; ++i)
	{
		if (!indeces[i].empty())
		{
			MPI_Isend(&indeces[i][0], static_cast<int>(indeces[i].size()), MPI_UINT32_T, static_cast<int>(i), 0, MPI_COMM_WORLD, &req);
			MPI_Isend(&ax[i][0], static_cast<int>(ax[i].size()), MPI_DOUBLE, static_cast<int>(i), 1, MPI_COMM_WORLD, &req);
			MPI_Isend(&ay[i][0], static_cast<int>(ay[i].size()), MPI_DOUBLE, static_cast<int>(i), 2, MPI_COMM_WORLD, &req);
			MPI_Isend(&az[i][0], static_cast<int>(az[i].size()), MPI_DOUBLE, static_cast<int>(i), 3, MPI_COMM_WORLD, &req);
		}
	}
	MPI_Status status;
	vector<vector<double > > recv_x(static_cast<size_t>(ntalkwithme)), recv_y(static_cast<size_t>(ntalkwithme)),
		recv_z(static_cast<size_t>(ntalkwithme));
	vector<vector<uint32_t> > recv_h(static_cast<size_t>(ntalkwithme));
	int count = 0;
	for (size_t i = 0; i < static_cast<size_t>(4*ntalkwithme); ++i)
	{
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		vector<int>::const_iterator it = binary_find(talkwithme.begin(), talkwithme.end(), status.MPI_SOURCE);
		assert(it != talkwithme.end());
		if (status.MPI_TAG == 0)
		{
			MPI_Get_count(&status, MPI_UINT32_T, &count);
			recv_h[static_cast<size_t>(it - talkwithme.begin())].resize(static_cast<size_t>(count));
			MPI_Recv(&recv_h[static_cast<size_t>(it - talkwithme.begin())][0], count, MPI_UINT32_T, status.MPI_SOURCE, status.MPI_TAG,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Get_count(&status, MPI_DOUBLE, &count);
			if (status.MPI_TAG == 1)
			{
				recv_x[static_cast<size_t>(it - talkwithme.begin())].resize(static_cast<size_t>(count));
				MPI_Recv(&recv_x[static_cast<size_t>(it - talkwithme.begin())][0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if (status.MPI_TAG == 2)
			{
				recv_y[static_cast<size_t>(it - talkwithme.begin())].resize(static_cast<size_t>(count));
				MPI_Recv(&recv_y[static_cast<size_t>(it - talkwithme.begin())][0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if (status.MPI_TAG == 3)
			{
				recv_z[static_cast<size_t>(it - talkwithme.begin())].resize(static_cast<size_t>(count));
				MPI_Recv(&recv_z[static_cast<size_t>(it - talkwithme.begin())][0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	acc.resize(Norg);
	size_t r_counter = 0;
	for (size_t i = 0; i < static_cast<size_t>(ntalkwithme); ++i)
	{
		size_t nrecv = recv_h[i].size();
		for (size_t j = 0; j < nrecv; j++)
		{
			++r_counter;
			acc[static_cast<size_t>(recv_h[i][j])].Set(recv_x[i][j], recv_y[i][j], recv_z[i][j]);
		}
	}
	assert(r_counter == Norg);
#else
struct reb_particle pt = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
for (size_t i = 0; i < Norg; ++i)
{
	Vector3D const& point = tess.GetMeshPoint(i);
	pt.x = point.x - CM.x;
	pt.y = point.y - CM.y;
	pt.z = point.z - CM.z;
	pt.vx = 0;
	pt.vy = 0;
	pt.vz = 0;
	pt.m = cells[i].density*tess.GetVolume(i);
	pt.hash = static_cast<uint32_t>(i);
	reb_add(r, pt);
}

// Start the integration
reb_integrate(r, r->dt);

// Get acc
for (size_t i = 0; i < Norg; ++i)
{
	acc[r->particles[i].hash].x = r->particles[i].ax;
	acc[r->particles[i].hash].y = r->particles[i].ay;
	acc[r->particles[i].hash].z = r->particles[i].az;
}
#endif
	// Cleanup
	reb_free_simulation(r);

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
			double distance = std::pow(abs(myCM - tess.GetCellCM(neigh[j])), -3.0);
			double new_distance = std::pow(smoothlength[i] * smoothlength[i] + smoothlength[neigh[j]] * smoothlength[neigh[j]], -1.5);
			acc[i] += otherM*(myCM - tess.GetCellCM(neigh[j])) *(distance-new_distance);
		}
	}
}
