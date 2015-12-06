#include "mpi_commands.hpp"
#include "../misc/utils.hpp"
#ifdef RICH_MPI

void MPI_exchange_data(const Tessellation& tess, vector<ComputationalCell>& cells)
{
	const vector<int>& correspondents = tess.GetDuplicatedProcs();
	const vector<vector<int> >& duplicated_points = tess.GetDuplicatedPoints();
	vector<vector<ComputationalCell> > incoming(correspondents.size());
	vector<MPI_Request> req(2*correspondents.size());
	vector<vector<double> > tempsend(correspondents.size());
	vector<double> temprecv;
	double temp = 0;
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		tempsend[i] = list_serialize(VectorValues(cells, duplicated_points[i]));
		int size = static_cast<int>(tempsend.size());
		if(size==0)
			MPI_Isend(&temp, 1, MPI_DOUBLE, correspondents[i], 0, MPI_COMM_WORLD, &req[2 * i]);
		else
			MPI_Isend(&tempsend[i][0], size, MPI_DOUBLE, correspondents[i], 0, MPI_COMM_WORLD, &req[2*i]);
	}
	const vector<vector<int> >& ghost_indices = tess.GetGhostIndeces();
	cells.resize(tess.GetTotalPointNumber());
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_DOUBLE,&count);
		temprecv.resize(static_cast<size_t>(count));
		MPI_Recv(&temprecv[0], count, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (count > 1)
		{
			vector<ComputationalCell> tempcells = list_unserialize(temprecv, cells[0]);
			size_t location = static_cast<size_t>(std::find(correspondents.begin(), correspondents.end(), status.MPI_SOURCE) -
				correspondents.begin());
			if (location >= correspondents.size())
				throw UniversalError("Bad location in mpi exchange");
			for (size_t j = 0; j < tempcells.size(); ++j)
				cells.at(ghost_indices.at(location).at(j)) = tempcells[j];
		}
	}
	MPI_Waitall(static_cast<int>(correspondents.size()), &req[0], MPI_STATUSES_IGNORE);
}

void exchange_cells(const Tessellation& tess, vector<ComputationalCell>& cells)
{
	const vector<int>& correspondents = tess.GetDuplicatedProcs();
	const vector<vector<int> >& duplicated_points = tess.GetDuplicatedPoints();
	vector<vector<ComputationalCell> > incoming(correspondents.size());
	vector<boost::mpi::request> requests;
	for (size_t i = 0; i < correspondents.size(); ++i)
	{
		requests.push_back(world.isend(correspondents[i], 0, VectorValues(cells, duplicated_points[i])));
		requests.push_back(world.irecv(correspondents[i], 0, incoming[i]));
	}
	const vector<vector<int> >& ghost_indices = tess.GetGhostIndeces();
	boost::mpi::wait_all(requests.begin(), requests.end());
	cells.resize(tess.GetTotalPointNumber());
	for (size_t i = 0; i < incoming.size(); ++i)
	{
		for (size_t j = 0; j < incoming.at(i).size(); ++j)
		{
			cells.at(ghost_indices.at(i).at(j)) = incoming.at(i).at(j);
		}
	}
}

MPI_Request SendVectorCells(vector<ComputationalCell> const& cells, int dest, int tag)
{
	MPI_Request req;
	vector<double> tosend = list_serialize(cells);
	int size = static_cast<int>(tosend.size());
	MPI_Isend(&tosend[0], size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &req);
	return req;
}
MPI_Request RecvVectorCells(vector<ComputationalCell> &cells, int dest, int tag);

MPI_Request SendVectorExtensives(vector<Extensive> const& extensives, int dest, int tag);
MPI_Request RecvVectorExtensives(vector<Extensive> &extensives, int dest, int tag);

MPI_Request SendVectorVector2D(vector<Vector2D> const& points, int dest, int tag);
MPI_Request RecvVectorVector2D(vector<Vector2D> &points, int dest, int tag);


#endif//RICH_MPI