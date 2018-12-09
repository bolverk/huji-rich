#include "hdf_write.hpp"
#include "../../misc/hdf5_utils.hpp"
#include "../../misc/int2str.hpp"

using namespace H5;

namespace
{

	template<class T> vector<T> read_vector_from_hdf5
	(const Group& file,
		const string& caption,
		const DataType& datatype)
	{
		DataSet dataset = file.openDataSet(caption);
		DataSpace filespace = dataset.getSpace();
		hsize_t dims_out[2];
		filespace.getSimpleExtentDims(dims_out, NULL);
		const size_t NX = static_cast<size_t>(dims_out[0]);
		vector<T> result(NX);
		dataset.read(&result[0], datatype);
		return result;
	}

	vector<double> read_double_vector_from_hdf5
	(Group& file, string const& caption)
	{
		return read_vector_from_hdf5<double>
			(file,
				caption,
				PredType::NATIVE_DOUBLE);
	}

	vector<int> read_int_vector_from_hdf5
	(const Group& file,
		const string& caption)
	{
		return read_vector_from_hdf5<int>
			(file,
				caption,
				PredType::NATIVE_INT);
	}

	vector<size_t> read_sizet_vector_from_hdf5
	(const Group& file,
		const string& caption)
	{
		return read_vector_from_hdf5<size_t>
			(file,
				caption,
				PredType::NATIVE_ULLONG);
	}
}

DiagnosticAppendix3D::~DiagnosticAppendix3D() {}

Snapshot3D::Snapshot3D(void) :
	mesh_points(),
#ifdef RICH_MPI
	proc_points(),
#endif
	volumes(),
	cells(),
	time(),
	cycle(),
	tracerstickernames() {}

Snapshot3D::Snapshot3D(const Snapshot3D& source) :
	mesh_points(source.mesh_points),
#ifdef RICH_MPI
	proc_points(source.proc_points),
#endif
	volumes(source.volumes),
	cells(source.cells),
	time(source.time),
	cycle(source.cycle),
	tracerstickernames(source.tracerstickernames)
{}

void WriteVoronoi(Voronoi3D const& tri, std::string const& filename)
{
	H5File file(H5std_string(filename), H5F_ACC_TRUNC);
	vector<double> x, y, z, vx, vy, vz;
	vector<size_t> Nfaces;
	vector<size_t> Nvert;
	vector<size_t> FacesInCell;
	vector<size_t> VerticesInFace;
	size_t Npoints = tri.GetPointNo();
	for (size_t i = 0; i < Npoints; ++i)
	{
		x.push_back(tri.GetMeshPoint(i).x);
		y.push_back(tri.GetMeshPoint(i).y);
		z.push_back(tri.GetMeshPoint(i).z);
	}
	write_std_vector_to_hdf5(file, x, "mesh_point_x");
	write_std_vector_to_hdf5(file, y, "mesh_point_y");
	write_std_vector_to_hdf5(file, z, "mesh_point_z");
	for (size_t i = 0; i < Npoints; ++i)
	{
		Nfaces.push_back(tri.GetCellFaces(i).size());
		for (size_t j = 0; j < Nfaces.back(); ++j)
			FacesInCell.push_back(tri.GetCellFaces(i)[j]);
	}
	IntType datatype(PredType::NATIVE_ULLONG);
	datatype.setOrder(H5T_ORDER_LE);
	write_std_vector_to_hdf5(file, Nfaces, "Number_of_faces_in_cell", datatype);
	write_std_vector_to_hdf5(file, FacesInCell, "Faces_in_cell", datatype);
	Npoints = tri.GetFacePoints().size();
	for (size_t i = 0; i < Npoints; ++i)
	{
		vx.push_back(tri.GetFacePoints()[i].x);
		vy.push_back(tri.GetFacePoints()[i].y);
		vz.push_back(tri.GetFacePoints()[i].z);
	}
	write_std_vector_to_hdf5(file, vx, "vertice_x");
	write_std_vector_to_hdf5(file, vy, "vertice_y");
	write_std_vector_to_hdf5(file, vz, "vertice_z");
	Npoints = tri.GetTotalFacesNumber();
	for (size_t i = 0; i < Npoints; ++i)
	{
		Nvert.push_back(tri.GetPointsInFace(i).size());
		for (size_t j = 0; j < Nvert.back(); ++j)
			VerticesInFace.push_back(tri.GetPointsInFace(i)[j]);
	}
	write_std_vector_to_hdf5(file, Nvert, "Number_of_vertices_in_face", datatype);
	write_std_vector_to_hdf5(file, VerticesInFace, "Vertices_in_face", datatype);
}

void WriteSnapshot3D(HDSim3D const& sim, std::string const& filename,const vector<DiagnosticAppendix3D*>& appendices)
{
	H5File file(H5std_string(filename), H5F_ACC_TRUNC);
	vector<ComputationalCell3D> const& cells = sim.getCells();
	Tessellation3D const& tess = sim.getTesselation();
	TracerStickerNames tsn = sim.GetTracerStickerNames();

	size_t Ncells = tess.GetPointNo();

	vector<double> temp(Ncells);

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tess.GetMeshPoint(i).x;
	write_std_vector_to_hdf5(file, temp, "X");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tess.GetMeshPoint(i).y;
	write_std_vector_to_hdf5(file, temp, "Y");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tess.GetMeshPoint(i).z;
	write_std_vector_to_hdf5(file, temp, "Z");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tess.GetCellCM(i).x;
	write_std_vector_to_hdf5(file, temp, "CMx");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tess.GetCellCM(i).y;
	write_std_vector_to_hdf5(file, temp, "CMy");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tess.GetCellCM(i).z;
	write_std_vector_to_hdf5(file, temp, "CMz");

#ifdef RICH_MPI
	// MPI
	Group mpi = file.createGroup("/mpi");
	Tessellation3D const& tproc = sim.getProcTesselation();
	Ncells = tproc.GetPointNo();
	temp.resize(Ncells);
	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tproc.GetMeshPoint(i).x;
	write_std_vector_to_hdf5(mpi, temp, "proc_X");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tproc.GetMeshPoint(i).y;
	write_std_vector_to_hdf5(mpi, temp, "proc_Y");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tproc.GetMeshPoint(i).z;
	write_std_vector_to_hdf5(mpi, temp, "proc_Z");
	Ncells = tess.GetPointNo();
	temp.resize(Ncells);
#endif

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].density;
	write_std_vector_to_hdf5(file, temp, "Density");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].pressure;
	write_std_vector_to_hdf5(file, temp, "Pressure");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].internal_energy;
	write_std_vector_to_hdf5(file, temp, "InternalEnergy");

	std::vector<size_t> ids(Ncells);
	for (size_t i = 0; i < Ncells; ++i)
		ids[i] = cells[i].ID;
	write_std_vector_to_hdf5(file, ids, "ID");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].velocity.x;
	write_std_vector_to_hdf5(file, temp, "Vx");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].velocity.y;
	write_std_vector_to_hdf5(file, temp, "Vy");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].velocity.z;
	write_std_vector_to_hdf5(file, temp, "Vz");

	Group tracers = file.createGroup("/tracers");
	Group stickers = file.createGroup("/stickers");

	for (size_t j = 0; j < tsn.tracer_names.size(); ++j)
	{
		for (size_t i = 0; i < Ncells; ++i)
			temp[i] = cells[i].tracers[j];
		write_std_vector_to_hdf5(tracers, temp, tsn.tracer_names[j]);
	}

	for (size_t j = 0; j <tsn.sticker_names.size(); ++j)
	{
		for (size_t i = 0; i < Ncells; ++i)
			temp[i] = cells[i].stickers[j];
		write_std_vector_to_hdf5(stickers, temp, tsn.sticker_names[j]);
	}

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = tess.GetVolume(i);
	write_std_vector_to_hdf5(file, temp, "Volume");

	vector<double> time(1, sim.GetTime());
	write_std_vector_to_hdf5(file, time, "Time");

	vector<int> cycle(1, static_cast<int>(sim.GetCycle()));
	write_std_vector_to_hdf5(file, cycle, "Cycle");

	// Appendices
	for (size_t i = 0; i < appendices.size(); ++i)
		write_std_vector_to_hdf5(file,(*(appendices.at(i)))(sim),appendices.at(i)->getName());
}

Snapshot3D ReadSnapshot3D
(const string& fname
#ifdef RICH_MPI
	, bool mpioverride
#endif
)
{
	Snapshot3D res;
	H5File file(fname, H5F_ACC_RDONLY);

	// Mesh points
	{
		const vector<double> x = read_double_vector_from_hdf5(file, "X");
		const vector<double> y = read_double_vector_from_hdf5(file, "Y");
		const vector<double> z = read_double_vector_from_hdf5(file, "Z");
		res.mesh_points.resize(x.size());
		for (size_t i = 0; i < x.size(); ++i)
			res.mesh_points.at(i) = Vector3D(x.at(i), y.at(i), z.at(i));
	}

#ifdef RICH_MPI
	// MPI
	if (!mpioverride)
	{
		Group mpi = file.openGroup("mpi");
		const vector<double> px = read_double_vector_from_hdf5(mpi, "proc_X");
		const vector<double> py = read_double_vector_from_hdf5(mpi, "proc_Y");
		const vector<double> pz = read_double_vector_from_hdf5(mpi, "proc_Z");
		res.proc_points.resize(px.size());
		for (size_t i = 0; i < px.size(); ++i)
			res.proc_points.at(i) = Vector3D(px.at(i), py.at(i), pz.at(i));
	}
#endif

	// Hydrodynamic
	{
		const vector<double> density = read_double_vector_from_hdf5(file, "Density");
		const vector<double> pressure = read_double_vector_from_hdf5(file, "Pressure");
		const vector<double> energy = read_double_vector_from_hdf5(file, "InternalEnergy");
		vector<size_t> IDs(density.size(), 0);
		ssize_t objcount =  file.getNumObjs();
		for (ssize_t i = 0; i < objcount; ++i)
		{
			std::string name = file.getObjnameByIdx(i);
			if (name.compare(std::string("ID"))==0)
				IDs = read_sizet_vector_from_hdf5(file, "ID");
		}
		const vector<double> x_velocity = read_double_vector_from_hdf5(file, "Vx");
		const vector<double> y_velocity = read_double_vector_from_hdf5(file, "Vy");
		const vector<double> z_velocity = read_double_vector_from_hdf5(file, "Vz");

		Group g_tracers = file.openGroup("tracers");
		Group g_stickers = file.openGroup("stickers");
		vector<vector<double> > tracers(g_tracers.getNumObjs());
		vector<string> tracernames(tracers.size());
		for (hsize_t n = 0; n < g_tracers.getNumObjs(); ++n)
		{
			const H5std_string name = g_tracers.getObjnameByIdx(n);
			tracernames[n] = name;
			tracers[n] = read_double_vector_from_hdf5(g_tracers, name);
		}

		vector<vector<int> > stickers(g_stickers.getNumObjs());
		vector<string> stickernames(stickers.size());
		for (hsize_t n = 0; n < g_stickers.getNumObjs(); ++n)
		{
			const H5std_string name = g_stickers.getObjnameByIdx(n);
			stickernames[n] = name;
			stickers[n] = read_int_vector_from_hdf5(g_stickers, name);
		}
		res.tracerstickernames.sticker_names = stickernames;
		res.tracerstickernames.tracer_names = tracernames;
		res.cells.resize(density.size());
		for (size_t i = 0; i < res.cells.size(); ++i)
		{
			res.cells.at(i).density = density.at(i);
			res.cells.at(i).pressure = pressure.at(i);
			res.cells.at(i).internal_energy = energy.at(i);
			res.cells.at(i).ID = IDs[i];
			res.cells.at(i).velocity.x = x_velocity.at(i);
			res.cells.at(i).velocity.y = y_velocity.at(i);
			res.cells.at(i).velocity.z = z_velocity.at(i);
			//res.cells.at(i).tracers.resize(tracernames.size());
			for (size_t j = 0; j < tracernames.size(); ++j)
				res.cells.at(i).tracers.at(j) = tracers.at(j).at(i);
			//res.cells.at(i).stickers.resize(stickernames.size());
			for (size_t j = 0; j < stickernames.size(); ++j)
				res.cells.at(i).stickers.at(j) = stickers.at(j).at(i) == 1;
		}
	}

	// Misc
	{
		const vector<double> time = read_double_vector_from_hdf5(file, "Time");
		res.time = time.at(0);

		const vector<double> volume = read_double_vector_from_hdf5(file, "Volume");
		res.volumes = volume;
		const vector<int> cycle = read_int_vector_from_hdf5(file, "Cycle");
		res.cycle = cycle.at(0);
	}

	return res;
}

#ifdef RICH_MPI
Snapshot3D ReDistributeData3D(string const& filename, Tessellation3D const& proctess, size_t snapshot_number)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// Read the data
	Snapshot3D snap;
	for (int i = 0; i < static_cast<int>(snapshot_number); ++i)
	{
		Snapshot3D temp = ReadSnapshot3D(filename + int2str(i) + ".h5", false);
		if (i == 0)
		{
			snap.time = temp.time;
			snap.cycle = temp.cycle;
			snap.tracerstickernames = temp.tracerstickernames;
		}
		size_t N = temp.cells.size();
		for (size_t i = 0; i < N; ++i)
		{
			if (PointInPoly(proctess, temp.mesh_points[i], static_cast<size_t>(rank)))
			{
				snap.cells.push_back(temp.cells[i]);
				snap.mesh_points.push_back(temp.mesh_points[i]);
			}
		}
	}
	return snap;
}
#endif