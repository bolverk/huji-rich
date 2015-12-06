#ifndef MPI_COMMANDS_HPP
#define MPI_COMMANDS_HPP 1

#ifdef RICH_MPI
#include <vector>
#include <mpi.h>
#include "../newtonian/two_dimensional/computational_cell_2d.hpp"
#include "../newtonian/two_dimensional/extensive.hpp"
#include "../tessellation/tessellation.hpp"
#include "../misc/serializable.hpp"
#include "../misc/utils.hpp"

using std::vector;

/*!
\brief Sends and revs data
\param tess The tessellation
\param cells The data to send/recv
\param ghost_or_sent True for ghost cells false for sent cells.
*/
template<class T>
void MPI_exchange_data(const Tessellation& tess, vector<T>& cells,bool ghost_or_sent)
{
	vector<int> correspondents;
	vector<vector<int> > duplicated_points;
	if (ghost_or_sent)
	{
		correspondents = tess.GetDuplicatedProcs();
		duplicated_points = tess.GetDuplicatedPoints();
	}
	else
	{
		correspondents = tess.GetSentProcs();
		duplicated_points = tess.GetSentPoints();
	}
	vector<MPI_Request> req(correspondents.size());
	vector<vector<double> > tempsend(correspondents.size());
	vector<double> temprecv;
	double temp = 0;
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		bool isempty = duplicated_points[i].empty();
		if(!isempty)
			tempsend[i] = list_serialize(VectorValues(cells, duplicated_points[i]));
		int size = static_cast<int>(tempsend[i].size());
		if (size == 0)
			MPI_Isend(&temp, 1, MPI_DOUBLE, correspondents[i], 1, MPI_COMM_WORLD, &req[i]);
		else
			MPI_Isend(&tempsend[i][0], size, MPI_DOUBLE, correspondents[i], 0, MPI_COMM_WORLD, &req[i]);
	}
	const vector<vector<int> >& ghost_indices = tess.GetGhostIndeces();
	if (ghost_or_sent)
		cells.resize(tess.GetTotalPointNumber());
	else
		cells = VectorValues(cells, tess.GetSelfPoint());
	vector<vector<T> > torecv(correspondents.size());
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		temprecv.resize(static_cast<size_t>(max(count,1)));
		MPI_Recv(&temprecv[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (status.MPI_TAG == 0)
		{
			size_t location = static_cast<size_t>(std::find(correspondents.begin(), correspondents.end(), status.MPI_SOURCE) -
				correspondents.begin());
			if (location >= correspondents.size())
				throw UniversalError("Bad location in mpi exchange");
			torecv[location] = list_unserialize(temprecv, cells[0]);
		}
	}
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		if (ghost_or_sent)
		{
			for (size_t j = 0; j < torecv[i].size(); ++j)
				cells.at(ghost_indices.at(i).at(j)) = torecv[i][j];
		}
		else
		{
			for (size_t j = 0; j < torecv[i].size(); ++j)
				cells.push_back(torecv[i][j]);
		}
	}
	MPI_Waitall(static_cast<int>(correspondents.size()), &req[0], MPI_STATUSES_IGNORE);
}

#endif //RICH_MPI
#endif // MPI_COMMANDS_HPP

