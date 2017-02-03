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
	if (cells.empty())
		throw UniversalError("Empty cell vector in MPI_exchange_data");
	T example_cell = cells[0];
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
			MPI_Isend(&temp, 1, MPI_DOUBLE, correspondents[i], 4, MPI_COMM_WORLD, &req[i]);
		else
			MPI_Isend(&tempsend[i][0], size, MPI_DOUBLE, correspondents[i], 5, MPI_COMM_WORLD, &req[i]);
	}
	const vector<vector<int> >& ghost_indices = tess.GetGhostIndeces();
	if (ghost_or_sent)
		cells.resize(tess.GetTotalPointNumber(),cells[0]);
	else
		cells = VectorValues(cells, tess.GetSelfPoint());
	vector<vector<T> > torecv(correspondents.size());
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		temprecv.resize(static_cast<size_t>(count));
		MPI_Recv(&temprecv[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (status.MPI_TAG == 5)
		{
			size_t location = static_cast<size_t>(std::find(correspondents.begin(), correspondents.end(), status.MPI_SOURCE) -
				correspondents.begin());
			if (location >= correspondents.size())
				throw UniversalError("Bad location in mpi exchange");
			torecv[location] = list_unserialize(temprecv, example_cell);
		}
		else
		{
			if (status.MPI_TAG != 4)
				throw UniversalError("Recv bad mpi tag");
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
	MPI_Barrier(MPI_COMM_WORLD);
}
/*!
\brief Sends and revs data
\param totalkwith The cpus to talk with
\param tosend The indeces in data to send ordered by cpu
\param cells The data to send
\return Th recv data ordered by cpu
*/
template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith,vector<vector<int> > const& tosend,
	vector<T>const& cells)
{
	vector<MPI_Request> req(totalkwith.size());
	vector<vector<double> > tempsend(totalkwith.size());
	vector<double> temprecv;
	double temp = 0;
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		bool isempty = tosend[i].empty();
		if (!isempty)
			tempsend[i] = list_serialize(VectorValues(cells, tosend[i]));
		int size = static_cast<int>(tempsend[i].size());
		if (size == 0)
			MPI_Isend(&temp, 1, MPI_DOUBLE, totalkwith[i], 4, MPI_COMM_WORLD, &req[i]);
		else
			MPI_Isend(&tempsend[i][0], size, MPI_DOUBLE, totalkwith[i], 5, MPI_COMM_WORLD, &req[i]);
	}
	vector<vector<T> > torecv(totalkwith.size());
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		temprecv.resize(static_cast<size_t>(count));
		MPI_Recv(&temprecv[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (status.MPI_TAG == 5)
		{
			size_t location = static_cast<size_t>(std::find(totalkwith.begin(), totalkwith.end(), status.MPI_SOURCE) -
				totalkwith.begin());
			if (location >= totalkwith.size())
				throw UniversalError("Bad location in mpi exchange");
			torecv[location] = list_unserialize(temprecv, cells[0]);
		}
		else
		{
			if (status.MPI_TAG != 4)
				throw UniversalError("Recv bad mpi tag");
		}
	}
	MPI_Waitall(static_cast<int>(totalkwith.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return torecv;
}

/*!
\brief Sends and revs data
\param totalkwith The cpus to talk with
\param cells The data to send ordered bys cpu
\param example An example cell for serialization
\return The recv data ordered by cpu
*/
template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<T> > const& cells,T const& example)
{
	vector<MPI_Request> req(totalkwith.size());
	vector<vector<double> > tempsend(totalkwith.size());
	vector<double> temprecv;
	double temp = 0;
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		bool isempty = cells[i].empty();
		if (!isempty)
			tempsend[i] = list_serialize(cells[i]);
		int size = static_cast<int>(tempsend[i].size());
		if (size == 0)
			MPI_Isend(&temp, 1, MPI_DOUBLE, totalkwith[i], 4, MPI_COMM_WORLD, &req[i]);
		else
			MPI_Isend(&tempsend[i][0], size, MPI_DOUBLE, totalkwith[i], 5, MPI_COMM_WORLD, &req[i]);
	}
	vector<vector<T> > torecv(totalkwith.size());
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		temprecv.resize(static_cast<size_t>(count));
		MPI_Recv(&temprecv[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (status.MPI_TAG == 5)
		{
			size_t location = static_cast<size_t>(std::find(totalkwith.begin(), totalkwith.end(), status.MPI_SOURCE) -
				totalkwith.begin());
			if (location >= totalkwith.size())
				throw UniversalError("Bad location in mpi exchange");
			torecv[location] = list_unserialize(temprecv, example);
		}
		else
		{
			if (status.MPI_TAG != 4)
				throw UniversalError("Recv bad mpi tag");
		}
	}
	MPI_Waitall(static_cast<int>(totalkwith.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return torecv;
}


/*!
\brief Sends and recvs data
\param tess The tessellation
\param cells The data to send/recv outermost vector defines the processors.
\param example An example type T, used in order to serialize.
\return The recv data ordered by cpu
*/
template<class T>
vector<vector<vector<T> > > MPI_exchange_data(const Tessellation& tess, vector<vector<vector<T> > > const& cells, T const& example)
{
	vector<int> correspondents = tess.GetDuplicatedProcs();
	vector<vector<vector<T> > > res(correspondents.size());
	vector<MPI_Request> req(2*correspondents.size());
	vector<vector<double> > tempsend(correspondents.size());
	vector<vector<int> > send_sizes(tempsend.size());
	vector<double> temprecv;
	vector<int> tempirecv;
	double temp = 0;
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		bool isempty = cells[i].empty();
		if (!isempty)
		{
			send_sizes[i].reserve(cells[i].size());
			tempsend[i].reserve(cells[i].size()*cells[i][0][0].getChunkSize());
			for (size_t j = 0; j < cells[i].size(); ++j)
			{
				send_sizes[i].push_back(static_cast<int>(cells[i][j].size()));
				vector<double> dtemp = list_serialize(cells[i][j]);
				tempsend[i].insert(tempsend[i].end(), dtemp.begin(), dtemp.end());
			}
		}
		else
			send_sizes[i].push_back(-1);
		int size = static_cast<int>(tempsend[i].size());
		if (size == 0)
		{
			MPI_Isend(&temp, 1, MPI_DOUBLE, correspondents[i], 4, MPI_COMM_WORLD, &req[2 * i]);
			MPI_Isend(&send_sizes[i][0], 1, MPI_INT, correspondents[i], 5, MPI_COMM_WORLD, &req[2 * i + 1]);
		}
		else
		{
			MPI_Isend(&tempsend[i][0], size, MPI_DOUBLE, correspondents[i], 6, MPI_COMM_WORLD, &req[2 * i]);
			MPI_Isend(&send_sizes[i][0], static_cast<int>(send_sizes[i].size()), MPI_INT, correspondents[i], 7, MPI_COMM_WORLD, &req[2 * i + 1]);
		}
	}
	vector<vector<double> > drecv(correspondents.size());
	vector<vector<int> > irecv(correspondents.size());
	for (size_t i = 0; i < 2*correspondents.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		if (status.MPI_TAG % 2 == 0)
		{
			MPI_Get_count(&status, MPI_DOUBLE, &count);
			size_t location = static_cast<size_t>(std::find(correspondents.begin(), correspondents.end(), status.MPI_SOURCE) -
				correspondents.begin());
			if (location >= correspondents.size())
				throw UniversalError("Bad location in mpi exchange");
			drecv[location].resize(static_cast<size_t>(count));
			MPI_Recv(&drecv[location][0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Get_count(&status, MPI_INT, &count);
			size_t location = static_cast<size_t>(std::find(correspondents.begin(), correspondents.end(), status.MPI_SOURCE) -
				correspondents.begin());
			if (location >= correspondents.size())
				throw UniversalError("Bad location in mpi exchange");
			irecv[location].resize(static_cast<size_t>(count));
			MPI_Recv(&irecv[location][0], count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		size_t counter = 0;
		for (size_t j = 0; j < irecv[i].size(); ++j)
		{
			if (irecv[i][j] > 0)
			{
				size_t to_read = static_cast<size_t>(irecv[i][j] * example.getChunkSize());
				vector<double> dtemp(drecv[i].begin() + counter, drecv[i].begin() + counter + to_read);
				counter += to_read;
				res[i].push_back(list_unserialize(dtemp, example));
			}
		}
	}
	MPI_Waitall(static_cast<int>(2*correspondents.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return res;
}

/*!
\brief Sends and recvs data
\param tess The tessellation
\param cells The data to send/recv outermost vector defines the processors.
\return The recv data ordered by cpu
*/
vector<vector<vector<int> > > MPI_exchange_data(const Tessellation& tess, vector<vector<vector<int> > > const& cells);

vector<vector<double> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<double> > &tosend);

vector<vector<int> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<int> > &tosend);

#endif //RICH_MPI
#endif // MPI_COMMANDS_HPP

