#include "hdf_write.hpp"
#include "hdf5_utils.hpp"

using namespace H5;

Snapshot::Snapshot(void) :
	mesh_points(){}

Snapshot::Snapshot(const Snapshot& source) :
	mesh_points(source.mesh_points){}

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
	write_std_vector_to_hdf5(file, Nfaces, "Number_of_faces_in_cell",datatype);
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

void WriteSnapshot(HDSim3D const& sim, std::string const& filename)
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
		temp[i] = cells[i].density;
	write_std_vector_to_hdf5(file, temp, "Density");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].pressure;
	write_std_vector_to_hdf5(file, temp, "Pressure");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].velocity.x;
	write_std_vector_to_hdf5(file, temp, "Vx");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].velocity.y;
	write_std_vector_to_hdf5(file, temp, "Vy");

	for (size_t i = 0; i < Ncells; ++i)
		temp[i] = cells[i].velocity.z;
	write_std_vector_to_hdf5(file, temp, "Vz");

	for (size_t j = 0; j < cells[0].tracers.size(); ++j)
	{
		for (size_t i = 0; i < Ncells; ++i)
			temp[i] = cells[i].tracers[j];
		write_std_vector_to_hdf5(file, temp, tsn.tracer_names[j]);
	}

	for (size_t j = 0; j < cells[0].stickers.size(); ++j)
	{
		for (size_t i = 0; i < Ncells; ++i)
			temp[i] = cells[i].stickers[j];
		write_std_vector_to_hdf5(file, temp, tsn.sticker_names[j]);
	}

	vector<double> time(1, sim.GetTime());
	write_std_vector_to_hdf5(file,time, "Time");

	vector<int> cycle(1, static_cast<int>(sim.GetCycle()));
	write_std_vector_to_hdf5(file, cycle, "Cycle");
}