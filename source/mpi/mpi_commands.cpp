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

#endif