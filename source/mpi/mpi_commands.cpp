#include "mpi_commands.hpp"

#ifdef RICH_MPI

vector<vector<double> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<double> > &tosend)
{
	vector<MPI_Request> req(totalkwith.size());
	vector<double> temprecv;
	double temp = 0;
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		int size = static_cast<int>(tosend[i].size());
		if (size == 0)
			MPI_Isend(&temp, 1, MPI_DOUBLE, totalkwith[i], 6, MPI_COMM_WORLD, &req[i]);
		else
			MPI_Isend(&tosend[i][0], size, MPI_DOUBLE, totalkwith[i], 7, MPI_COMM_WORLD, &req[i]);
	}
	vector<vector<double> > torecv(totalkwith.size());
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_DOUBLE, &count);
		temprecv.resize(static_cast<size_t>(count));
		MPI_Recv(&temprecv[0], count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (status.MPI_TAG == 7)
		{
			size_t location = static_cast<size_t>(std::find(totalkwith.begin(), totalkwith.end(), status.MPI_SOURCE) -
				totalkwith.begin());
			if (location >= totalkwith.size())
				throw UniversalError("Bad location in mpi exchange");
			torecv[location] = temprecv;
		}
		else
		{
			if (status.MPI_TAG != 6)
				throw UniversalError("Recv bad mpi tag");
		}
	}
	MPI_Waitall(static_cast<int>(totalkwith.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return torecv;
}


vector<vector<int> > MPI_exchange_data(const vector<int>& totalkwith, vector<vector<int> > &tosend)
{
	vector<MPI_Request> req(totalkwith.size());
	vector<int> temprecv;
	int temp = 0;
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		int size = static_cast<int>(tosend[i].size());
		if (size == 0)
			MPI_Isend(&temp, 1, MPI_INT, totalkwith[i], 8, MPI_COMM_WORLD, &req[i]);
		else
			MPI_Isend(&tosend[i][0], size, MPI_INT, totalkwith[i], 9, MPI_COMM_WORLD, &req[i]);
	}
	vector<vector<int> > torecv(totalkwith.size());
	for (size_t i = 0; i < totalkwith.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_INT, &count);
		temprecv.resize(static_cast<size_t>(count));
		MPI_Recv(&temprecv[0], count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (status.MPI_TAG == 9)
		{
			size_t location = static_cast<size_t>(std::find(totalkwith.begin(), totalkwith.end(), status.MPI_SOURCE) -
				totalkwith.begin());
			if (location >= totalkwith.size())
				throw UniversalError("Bad location in mpi exchange");
			torecv[location] = temprecv;
		}
		else
		{
			if (status.MPI_TAG != 8)
				throw UniversalError("Recv bad mpi tag");
		}
	}
	MPI_Waitall(static_cast<int>(totalkwith.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return torecv;
}

void MPI_exchange_data2(const Tessellation& tess, vector<double>& cells, bool ghost_or_sent)
{
	if (cells.empty())
		throw UniversalError("Empty double cell vector in MPI_exchange_data");
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
	vector<vector<double> > tempsend(correspondents.size());
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		bool isempty = duplicated_points[i].empty();
		if (!isempty)
			tempsend[i] = VectorValues(cells, duplicated_points[i]);
	}
	vector<vector<double> > torecv = MPI_exchange_data(correspondents, tempsend);

	const vector<vector<int> >& ghost_indices = tess.GetGhostIndeces();
	if (ghost_or_sent)
		cells.resize(tess.GetTotalPointNumber(), cells[0]);
	else
		cells = VectorValues(cells, tess.GetSelfPoint());
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
}

vector<vector<vector<int> > > MPI_exchange_data(const Tessellation& tess, vector<vector<vector<int> > > const& cells)
{
	vector<int> correspondents = tess.GetDuplicatedProcs();
	vector<vector<vector<int> > > res(correspondents.size());
	vector<vector<int> > tempsend(correspondents.size());
	vector<MPI_Request> req(2 * correspondents.size());


	vector<vector<int> > send_sizes(correspondents.size());
	vector<int> temprecv;
	vector<int> tempsrecv;
	int temp = 0;
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		bool isempty = cells[i].empty();
		if (!isempty)
		{
			send_sizes[i].reserve(cells[i].size());
			for (size_t j = 0; j < cells[i].size(); ++j)
			{
				send_sizes[i].push_back(static_cast<int>(cells[i][j].size()));
				tempsend[i].insert(tempsend[i].end(), cells[i][j].begin(), cells[i][j].end());
			}
		}
		else
			send_sizes[i].push_back(-1);
		int size = static_cast<int>(tempsend[i].size());
		if (size == 0)
		{
			MPI_Isend(&temp, 1, MPI_INT, correspondents[i], 4, MPI_COMM_WORLD, &req[2 * i]);
			MPI_Isend(&send_sizes[i][0], 1, MPI_INT, correspondents[i], 5, MPI_COMM_WORLD, &req[2 * i + 1]);
		}
		else
		{
			MPI_Isend(&tempsend[i][0], size, MPI_INT, correspondents[i], 6, MPI_COMM_WORLD, &req[2 * i]);
			MPI_Isend(&send_sizes[i][0], static_cast<int>(send_sizes[i].size()), MPI_INT, correspondents[i], 7, MPI_COMM_WORLD, &req[2 * i + 1]);
		}
	}
	vector<vector<int> > srecv(correspondents.size());
	vector<vector<int> > irecv(correspondents.size());
	for (size_t i = 0; i < 2 * correspondents.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int count;
		if (status.MPI_TAG % 2 == 0)
		{
			MPI_Get_count(&status, MPI_INT, &count);
			size_t location = static_cast<size_t>(std::find(correspondents.begin(), correspondents.end(), status.MPI_SOURCE) -
				correspondents.begin());
			if (location >= correspondents.size())
				throw UniversalError("Bad location in mpi exchange");
			irecv[location].resize(static_cast<size_t>(count));
			MPI_Recv(&irecv[location][0], count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Get_count(&status, MPI_INT, &count);
			size_t location = static_cast<size_t>(std::find(correspondents.begin(), correspondents.end(), status.MPI_SOURCE) -
				correspondents.begin());
			if (location >= correspondents.size())
				throw UniversalError("Bad location in mpi exchange");
			srecv[location].resize(static_cast<size_t>(count));
			MPI_Recv(&srecv[location][0], count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		size_t counter = 0;
		for (size_t j = 0; j < srecv[i].size(); ++j)
		{
			if (srecv[i][j] < 1)
				continue;
			size_t size_add = static_cast<size_t>(srecv[i][j]);
			temprecv.assign(irecv[i].begin() + counter, irecv[i].begin() + counter + size_add);
			res[i].push_back(temprecv);
			counter += size_add;
		}
	}
	MPI_Waitall(static_cast<int>(2 * correspondents.size()), &req[0], MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	return res;
}

#endif
