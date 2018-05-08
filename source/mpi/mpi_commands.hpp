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
#include "stdint.h"
#include <boost/align.hpp>

using std::vector;

/*!
\brief Sends and revs data
\param tess The tessellation
\param cells The data to send/recv
\param ghost_or_sent True for ghost cells false for sent cells.
*/
template<class T>
void MPI_exchange_data(const Tessellation& tess, vector<T>& cells, bool ghost_or_sent);

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
	if (!req.empty())
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
void MPI_exchange_data(const Tessellation3D& tess, vector<T>& cells, bool ghost_or_sent);

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
	if (!req.empty())
		MPI_Waitall(static_cast<int>(correspondents.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
}

/*!
\brief Sends and revs data
\param totalkwith The cpus to talk with
\param tosend The indeces in data to send ordered by cpu
\param cells The data to send
\return The recv data ordered by cpu
*/
template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<int> > const& tosend,
	vector<T>const& cells);

template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith,vector<vector<int> > const& tosend,
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
	if (!req.empty())
		MPI_Waitall(static_cast<int>(totalkwith.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return torecv;
}


template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<int> > const& tosend,
	vector<T, boost::alignment::aligned_allocator<T, 32> >const& cells);
template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<int> > const& tosend,
	vector<T,boost::alignment::aligned_allocator<T,32> >const& cells)
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
	if (!req.empty())
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
	vector<T>const& cells);

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
	if(!req.empty())
		MPI_Waitall(static_cast<int>(totalkwith.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return torecv;
}

template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<size_t> > const& tosend,
	vector<T,boost::alignment::aligned_allocator<T,32> >const& cells);

template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<size_t> > const& tosend,
	vector<T, boost::alignment::aligned_allocator<T, 32> >const& cells)
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
	if (!req.empty())
		MPI_Waitall(static_cast<int>(totalkwith.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return torecv;
}



template <class T>
vector<vector<T> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<T> > const& tosend, T const& demo);

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
	if (!req.empty())
		MPI_Waitall(static_cast<int>(totalkwith.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return torecv;
}

template <class T>
vector<vector<vector<T> > > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<vector<T > > > const& tosend,
	T const& demo);

template <class T>
vector<vector<vector<T> > > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<vector<T > > > const& tosend,
	T const& demo)
{
	vector<vector<vector<T > > > res(totalkwith.size());
	vector<MPI_Request> req(2*totalkwith.size());
	vector<vector<double> > tempsend(2*totalkwith.size());
	vector<vector<int> > send_sizes(totalkwith.size());
	vector<T> temprecv;
	double temp = 0;
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		bool isempty = tosend[i].empty();
		if (!isempty)
		{
			send_sizes[i].reserve(tosend[i].size());
			for (size_t j = 0; j < tosend[i].size(); ++j)
			{
				send_sizes[i].push_back(static_cast<int>(tosend[i][j].size()));
				vector<double> dtemp = list_serialize(tosend[i][j]);
				tempsend[i].insert(tempsend[i].end(), dtemp.begin(), dtemp.end());
			}		
		}
		else
			send_sizes[i].push_back(-1);
		int size = static_cast<int>(tempsend[i].size());
		if (size == 0)
		{
			MPI_Isend(&temp, 1, MPI_DOUBLE, totalkwith[i], 4, MPI_COMM_WORLD, &req[2 * i]);
			MPI_Isend(&send_sizes[i][0], 1, MPI_INT, totalkwith[i], 5, MPI_COMM_WORLD, &req[2 * i + 1]);
		}
		else
		{
			MPI_Isend(&tempsend[i][0], size, MPI_DOUBLE, totalkwith[i], 6, MPI_COMM_WORLD, &req[2 * i]);
			MPI_Isend(&send_sizes[i][0], static_cast<int>(send_sizes[i].size()), MPI_INT, totalkwith[i], 7, MPI_COMM_WORLD, &req[2 * i + 1]);
		}
	}
	vector<vector<double> > srecv(totalkwith.size());
	vector<vector<int> > irecv(totalkwith.size());
	for (size_t i = 0; i < 2 * totalkwith.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		if (status.MPI_TAG % 2 == 1)
		{
			MPI_Get_count(&status, MPI_INT, &count);
			size_t location = static_cast<size_t>(std::find(totalkwith.begin(), totalkwith.end(), status.MPI_SOURCE) -
				totalkwith.begin());
			if (location >= totalkwith.size())
				throw UniversalError("Bad location in mpi exchange");
			irecv[location].resize(static_cast<size_t>(count));
			MPI_Recv(&irecv[location][0], count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Get_count(&status, MPI_DOUBLE, &count);
			size_t location = static_cast<size_t>(std::find(totalkwith.begin(), totalkwith.end(), status.MPI_SOURCE) -
				totalkwith.begin());
			if (location >= totalkwith.size())
				throw UniversalError("Bad location in mpi exchange");
			srecv[location].resize(static_cast<size_t>(count));
			MPI_Recv(&srecv[location][0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		if (irecv[i].at(0) < 1)
			continue;
		vector<T> T_temp = list_unserialize(srecv[i], demo);
		if (T_temp.empty())
			continue;
		size_t counter = 0;
		for (size_t j = 0; j < irecv[i].size(); ++j)
		{
			if (irecv[i][j] < 1)
				continue;
			size_t size_add = static_cast<size_t>(irecv[i][j]);
			vector<T> T_add(T_temp.begin() + counter, T_temp.begin() + counter + size_add);
			res[i].push_back(T_add);
			counter += size_add;
		}
	}
	if (!req.empty())
		MPI_Waitall(static_cast<int>(2*totalkwith.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return res;
}

void MPI_exchange_data2(const Tessellation3D& tess, vector<double>& cells, bool ghost_or_sent);

vector<vector<double> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<double> > &tosend);

vector<vector<int> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<int> > &tosend);

vector<vector<vector<int> > > MPI_exchange_data(const Tessellation3D& tess, vector<vector<vector<int> > > &tosend);

vector<vector<vector<size_t> > > MPI_exchange_data(const Tessellation3D& tess, vector<vector<vector<size_t> > > &tosend);

vector<vector<vector<double> > > MPI_exchange_data(const Tessellation3D& tess, vector<vector<vector<double> > > &tosend);

vector<vector<size_t> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<size_t> > &tosend);

void MPI_exchange_data(const Tessellation3D& tess, vector<char>& cells, bool ghost_or_sent);
#endif //RICH_MPI
#endif // MPI_COMMANDS_HPP

