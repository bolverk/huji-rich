#include "HilbertProcPositions.hpp"
#include "../3D/GeometryCommon/HilbertOrder3D.hpp"
#include "../misc/utils.hpp"
#include <limits>
#ifdef RICH_MPI
#include <mpi.h>
#endif

namespace
{
#ifdef RICH_MPI
	vector<unsigned long long> SortAllPoints(size_t Ncor, size_t Nproc, vector<size_t> &Hindeces, Tessellation3D const& tess)
	{
		assert(Hindeces.size() == Ncor && Ncor>0);
		double segfraction = 0.04;
		int Ncor2 = static_cast<int>(Ncor);
		int Ntotal;
		MPI_Allreduce(&Ncor2, &Ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		Ntotal = Ncor2;
		double Npart = static_cast<double>(Ntotal)*segfraction / static_cast<double>(Nproc);
		vector<Vector3D> cor(Ncor);
		for (size_t i = 0; i < Ncor; ++i)
			cor[i] = tess.GetMeshPoint(i);
		vector<size_t> sort_indeces(Ncor, 0);
		sort_index(Hindeces, sort_indeces);
		std::sort(Hindeces.begin(), Hindeces.end());
		vector<unsigned long long> segment_indeces;
		segment_indeces.push_back(static_cast<unsigned long long>(Hindeces[0]));
		size_t counter = 1;

		for (size_t i = 1; i < Ncor; ++i)
		{
			if (counter >= static_cast<size_t>(Npart))
			{
				segment_indeces.push_back(static_cast<unsigned long long>(Hindeces[i]));
				counter = 1;
			}
			else
				++counter;
		}
		if (counter > static_cast<size_t>(Npart*0.3))
		{
			segment_indeces.push_back(static_cast<unsigned long long>(Hindeces.back()));
		}
		else
		{
			segment_indeces.back() = Hindeces.back();
		}

		vector<int> nsegments_per_proc(Nproc), no_disp(Nproc, 0);
		int nsegments_local = static_cast<int>(segment_indeces.size());
		MPI_Allgather(&nsegments_local, 1, MPI_INT, &nsegments_per_proc[0], 1, MPI_INT, MPI_COMM_WORLD);
		int seg_size_total = 0;
		for (size_t i = 0; i < Nproc; ++i)
		{
			seg_size_total += nsegments_per_proc[i];
			if (i == 0)
				no_disp[i] = 0;
			else
				no_disp[i] = nsegments_per_proc[i - 1] + no_disp[i - 1];
		}

		vector<unsigned long long> all_segments(seg_size_total);
		MPI_Allgatherv(&segment_indeces[0], nsegments_local, MPI_UNSIGNED_LONG_LONG, &all_segments[0],
			&nsegments_per_proc[0], &no_disp[0], MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

		std::sort(all_segments.begin(), all_segments.end());
		vector<unsigned long long> res_segments(Nproc);
		double seg_per_proc = static_cast<double>(seg_size_total - 1) / static_cast<double>(Nproc);
		for (size_t i = 0; i < Nproc; ++i)
			res_segments[i] = all_segments.at(static_cast<size_t>(static_cast<double>(i + 1)*seg_per_proc));
		res_segments.back() = all_segments.back();

		return res_segments;
	}
	
	vector<Vector3D> GetNewProcPoints(vector<unsigned long long> const& Hxcor, Tessellation3D const& tess,
		vector<size_t> const& H_indeces, size_t Nproc)
	{
		vector<Vector3D> maxv(Nproc, Vector3D(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
			-std::numeric_limits<double>::max())), minv(Nproc, Vector3D(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
				std::numeric_limits<double>::max()));
		size_t Ncor = tess.GetPointNo();
		for (size_t i = 0; i < Ncor; ++i)
		{
			size_t index = static_cast<size_t>(std::lower_bound(Hxcor.begin(), Hxcor.end(), static_cast<unsigned long long>(H_indeces[i]))
				- Hxcor.begin());
			Vector3D const& cor = tess.GetMeshPoint(i);
			minv[index].x = std::min(minv[index].x, cor.x);
			minv[index].y = std::min(minv[index].y, cor.y);
			minv[index].z = std::min(minv[index].z, cor.z);
			maxv[index].x = std::max(maxv[index].x, cor.x);
			maxv[index].y = std::max(maxv[index].y, cor.y);
			maxv[index].z = std::max(maxv[index].z, cor.z);
		}

		vector<double> sendtemp(Nproc * 3, 0), recvtemp(Nproc * 3, 0);
		for (size_t i = 0; i < Nproc; ++i)
		{
			sendtemp[3 * i] = maxv[i].x;
			sendtemp[3 * i + 1] = maxv[i].y;
			sendtemp[3 * i + 2] = maxv[i].z;
		}
		MPI_Allreduce(&sendtemp[0], &recvtemp[0], static_cast<int>(3 * Nproc), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		for (size_t i = 0; i < Nproc; ++i)
		{
			maxv[i].x = recvtemp[3 * i];
			maxv[i].y = recvtemp[3 * i + 1];
			maxv[i].z = recvtemp[3 * i + 2];
		}
		for (size_t i = 0; i < Nproc; ++i)
		{
			sendtemp[3 * i] = minv[i].x;
			sendtemp[3 * i + 1] = minv[i].y;
			sendtemp[3 * i + 2] = minv[i].z;
		}

		MPI_Allreduce(&sendtemp[0], &recvtemp[0], static_cast<int>(3 * Nproc), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		for (size_t i = 0; i < Nproc; ++i)
		{
			minv[i].x = recvtemp[3 * i];
			minv[i].y = recvtemp[3 * i + 1];
			minv[i].z = recvtemp[3 * i + 2];
		}

		for (size_t i = 0; i < Nproc; ++i)
			maxv[i] = 0.5*(maxv[i] + minv[i]);

		return maxv;
	}
#endif

}

vector<Vector3D> HilbertProcPositions(Tessellation3D const & tess)
{
	size_t Ncor = tess.GetPointNo();
	vector<Vector3D> res(Ncor);
#ifdef RICH_MPI
	int ws = 32;
	MPI_Comm_size(MPI_COMM_WORLD, &ws);
	size_t Nproc = static_cast<size_t>(ws);
	int Ntotal = static_cast<int>(Ncor);
	int Ncor2 = static_cast<int>(Ncor);
	MPI_Allreduce(&Ncor2, &Ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	vector<Vector3D> cor(Ncor);
	for (size_t i = 0; i < Ncor; ++i)
		cor[i] = tess.GetMeshPoint(i);
	size_t Hmax;
	vector<size_t> H_indeces = GetGlobalHibertIndeces(cor, tess.GetBoxCoordinates().first, tess.GetBoxCoordinates().second, Hmax);
	vector<size_t> H_indeces2(H_indeces);
	vector<unsigned long long> Hindex = SortAllPoints(Ncor, Nproc, H_indeces, tess);
	res = GetNewProcPoints(Hindex, tess, H_indeces2, Nproc);
#endif
	return res;
}
