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
TreeSelfGravity::TreeSelfGravity(double smoothlength, size_t nx, size_t ny, size_t nz):smoothlength_(smoothlength),
nx_(nx),ny_(ny),nz_(nz){}

void TreeSelfGravity::operator()(const Tessellation3D & tess, const vector<ComputationalCell3D>& cells, 
	const vector<Conserved3D>& /*fluxes*/, const double /*time*/, TracerStickerNames const & /*tracerstickernames*/,
	vector<Vector3D>& acc) const
{
	std::pair<Vector3D,Vector3D> box = tess.GetBoxCoordinates();
	Vector3D CM = 0.5*(box.first + box.second);
	Vector3D diff = box.second - box.first;
	double boxsize = std::max(std::max(diff.x, diff.y), diff.z);

	struct reb_simulation* const r = reb_create_simulation();
	// Setup constants
	r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
	r->gravity = reb_simulation::REB_GRAVITY_TREE;
	r->boundary = reb_simulation::REB_BOUNDARY_NONE;
	r->opening_angle2 = 0.5;          // This constant determines the accuracy of the tree code gravity estimate.
	r->G = 1;
	r->softening = smoothlength_;         // Gravitational softening length
	r->dt = DBL_MIN*100;         // Timestep
	// Setup root boxes for gravity tree.
	// Here, we use 2x2=4 root boxes (each with length 'boxsize')
	// This allows you to use up to 4 MPI nodes.
	reb_configure_box(r, boxsize, static_cast<int>(nx_),static_cast<int>( ny_), static_cast<int>(nz_));

	// Get total number of points
	vector<int> point_num(1, static_cast<int> (tess.GetPointNo()));
	int Ntotal = point_num[0];
#ifdef RICH_MPI
	int ws = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	point_num.resize(static_cast<size_t>(ws));
	MPI_Gather(&Ntotal, 1, MPI_INT, &point_num[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	for (int i = 1; i < ws; ++i)
		Ntotal += point_num[static_cast<size_t>(i)];
#endif
	// Get all of the points and their masses
	vector<Vector3D> allpoints(static_cast<size_t>(Ntotal));
	vector<double> masses(allpoints.size());
	for (size_t i = 0; i < static_cast<size_t>(point_num[0]); ++i)
	{
		allpoints[i] = tess.GetMeshPoint(i);
		masses[i] = tess.GetVolume(i)*cells[i].density;
	}
#ifdef RICH_MPI
	vector<double> temp;
	vector<size_t> cum_num(point_num.size(),0);
	MPI_Request req,req2;
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank != 0)
	{
		temp = list_serialize(allpoints);
		MPI_Isend(&temp[0], static_cast<int>(temp.size()), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD,&req);
		MPI_Isend(&masses[0], static_cast<int>(masses.size()), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &req2);
	}
	else
	{
		vector<Vector3D> vtemp;
		for (size_t i = 1; i < static_cast<size_t>(ws); ++i)
			cum_num[i] = point_num[i - 1];
		for (size_t i = 1; i < static_cast<size_t>(ws); ++i)
			cum_num[i] += cum_num[i - 1];
		MPI_Status status;
		for (int i = 0; i < (2*ws-2); ++i)
		{
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			int count = 0;
			MPI_Get_count(&status, MPI_DOUBLE, &count);
			temp.resize(static_cast<size_t>(count));
			MPI_Recv(&temp[0], count, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (status.MPI_TAG == 1)
			{
				vtemp = list_unserialize(temp, tess.GetMeshPoint(0));
				std::copy(vtemp.begin(), vtemp.end(), allpoints.begin() +
					static_cast<size_t>(cum_num[static_cast<size_t>(status.MPI_SOURCE)]));
			}
			else
			{
				assert(status.MPI_TAG == 2);
				std::copy(temp.begin(), temp.end(), masses.begin() +
					static_cast<size_t>(cum_num[static_cast<size_t>(status.MPI_SOURCE)]));
			}
		}
	}

#endif

	size_t Npart = static_cast<size_t>(Ntotal);
	acc.resize(Npart);
#ifdef RICH_MPI
	if (rank == 0)
	{
#endif
	struct reb_particle pt = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };	
		for (size_t i = 0; i < Npart; ++i)
		{
			Vector3D const& point = allpoints[i];
			pt.x = point.x;
			pt.y = point.y;
			pt.z = point.z;
			pt.vx = 0;
			pt.vy = 0;
			pt.vz = 0;
			pt.m = masses[i];
			pt.hash = static_cast<uint32_t>(i);
			reb_add(r, pt);
		}

		// Start the integration
		reb_integrate(r, r->dt);

		// Get acc
		for (size_t i = 0; i < Npart; ++i)
		{
			acc[r->particles[i].hash].x = r->particles[i].ax;
			acc[r->particles[i].hash].y = r->particles[i].ay;
			acc[r->particles[i].hash].z = r->particles[i].az;
		}

#ifdef RICH_MPI
		// send acc back to cpus
		temp = list_serialize(acc);
		for (int i = 1; i < ws; ++i)
			MPI_Isend(&temp[3 * cum_num[static_cast<size_t>(i)]], 3 * point_num[static_cast<size_t>(i)], MPI_DOUBLE, i, 0,
				MPI_COMM_WORLD, &req);
	}
	else
	{
		// recv the acc
		temp.resize(3 * Npart);
		MPI_Recv(&temp[0], static_cast<int>(3 * Npart), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
		acc.resize(Npart);
#endif


	

	// Cleanup
	reb_free_simulation(r);
}
