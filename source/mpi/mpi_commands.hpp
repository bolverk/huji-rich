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
#include "../3D/GeometryCommon/Tessellation3D.hpp"

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
\param tess The tessellation
\param cells The data to send/recv
\param ghost_or_sent True for ghost cells false for sent cells.
*/
template<class T>
void MPI_exchange_data(const Tessellation3D& tess, vector<T>& cells, bool ghost_or_sent)
{
	if (cells.empty())
		throw UniversalError("Empty cell vector in MPI_exchange_data");
	T example_cell = cells[0];
	vector<int> correspondents;
	vector<vector<size_t> > duplicated_points;
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
		if (!isempty)
			tempsend[i] = list_serialize(VectorValues(cells, duplicated_points[i]));
		int size = static_cast<int>(tempsend[i].size());
		if (size == 0)
			MPI_Isend(&temp, 1, MPI_DOUBLE, correspondents[i], 4, MPI_COMM_WORLD, &req[i]);
		else
			MPI_Isend(&tempsend[i][0], size, MPI_DOUBLE, correspondents[i], 5, MPI_COMM_WORLD, &req[i]);
	}
	const vector<vector<size_t> >& ghost_indices = tess.GetGhostIndeces();
	if (ghost_or_sent)
		cells.resize(tess.GetTotalPointNumber(), cells[0]);
	else
		cells = VectorValues(cells, tess.GetSelfIndex());
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
\param tosend The indeces in data to send ordered by cpu
\param cells The data to send
\return The recv data ordered by cpu
*/
template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<size_t> > const& tosend,
	vector<T>const& cells)
{
	assert(!cells.empty());
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

template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<T> > const& tosend,T const& demo)
{
	vector<MPI_Request> req(totalkwith.size());
	vector<vector<double> > tempsend(totalkwith.size());
	vector<double> temprecv;
	double temp = 0;
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		bool isempty = tosend[i].empty();
		if (!isempty)
			tempsend[i] = list_serialize(tosend[i]);
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
			torecv[location] = list_unserialize(temprecv, demo);
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

vector<vector<double> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<double> > &tosend);

vector<vector<int> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<int> > &tosend);

vector<vector<vector<int> > > MPI_exchange_data(const Tessellation3D& tess, vector<vector<vector<int> > > const& tosend);

#endif //RICH_MPI
#endif // MPI_COMMANDS_HPP

