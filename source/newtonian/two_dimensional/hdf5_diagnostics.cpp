#include "hdf5_diagnostics.hpp"
#include "../../misc/hdf5_utils.hpp"

using namespace H5;

namespace {
  vector<double> read_double_vector_from_hdf5(H5File& file,string const& data_name)
  {
    DataSet dataset = file.openDataSet(data_name);
    DataSpace filespace = dataset.getSpace();
    hsize_t dims_out[2];
    filespace.getSimpleExtentDims(dims_out,NULL);
    int NX = static_cast<int>(dims_out[0]);
    vector<double> result(static_cast<size_t>(NX));
    dataset.read(&result[0],PredType::NATIVE_DOUBLE);
    return result;
  }

  vector<size_t> read_sizet_vector_from_hdf5(H5File& file,string const& data_name)
  {
    DataSet dataset = file.openDataSet(data_name);
    DataSpace filespace = dataset.getSpace();
    hsize_t dims_out[2];
    filespace.getSimpleExtentDims(dims_out,NULL);
    int NX = static_cast<int>(dims_out[0]);
    vector<unsigned> result(static_cast<size_t>(NX));
    dataset.read(&result[0],PredType::NATIVE_UINT);
    return list_static_cast<size_t,unsigned>(result);
  }

  vector<int> read_int_vector_from_hdf5(H5File& file,string const& data_name)
  {
    DataSet dataset = file.openDataSet(data_name);
    DataSpace filespace = dataset.getSpace();
    hsize_t dims_out[2];
    filespace.getSimpleExtentDims(dims_out,NULL);
    int NX=static_cast<int>(dims_out[0]);
    vector<int> result(static_cast<size_t>(NX));
    dataset.read(&result[0],PredType::NATIVE_INT);
    return result;
  }
}

namespace {
  class MeshGeneratingPointCoordinate: public Index2Member<double>
  {
  public:

    MeshGeneratingPointCoordinate(const Tessellation& tess,
				  double Vector2D::* component):
      tess_(tess), component_(component) {}

    size_t getLength(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
    }

    double operator()(size_t i) const
    {
      return tess_.GetMeshPoint(static_cast<int>(i)).*component_;
    }

  private:
    const Tessellation& tess_;
    double Vector2D::* component_;
  };

  class SingleCellPropertyExtractor
  {
  public:

    virtual double operator()(const Primitive& p) const = 0;

    virtual ~SingleCellPropertyExtractor(void) {}
  };

  class ThermalPropertyExtractor: public SingleCellPropertyExtractor
  {
  public:

    ThermalPropertyExtractor(double Primitive::* var):
      var_(var) {}

    double operator()(const Primitive& p) const
    {
      return p.*var_;
    }

  private:
    double Primitive::* var_;
  };

  class CellVelocityComponentExtractor: public SingleCellPropertyExtractor
  {
  public:

    CellVelocityComponentExtractor(double Vector2D::* component):
      component_(component) {}

    double operator()(const Primitive& p) const
    {
      return p.Velocity.*component_;
    }

  private:
    double Vector2D::* component_;
  };

  class CellsPropertyExtractor: public Index2Member<double>
  {
  public:

    CellsPropertyExtractor(const hdsim& sim,
			   const SingleCellPropertyExtractor& scpe):
      sim_(sim), scpe_(scpe) {}

    size_t getLength(void) const
    {
      return static_cast<size_t>(sim_.GetCellNo());
    }

    double operator()(size_t i) const
    {
      return scpe_(sim_.GetCell(static_cast<int>(i)));
    }

  private:
    const hdsim& sim_;
    const SingleCellPropertyExtractor& scpe_;
  };

  class ConvexHullData
  {
  public:
    
    vector<double> xvert;
    vector<double> yvert;
    vector<double> nvert;

    ConvexHullData(const Tessellation& tess):
      xvert(),
      yvert(),
      nvert(static_cast<size_t>(tess.GetPointNo()))
    {
      xvert.reserve(7*static_cast<size_t>(tess.GetPointNo()));
      yvert.reserve(7*static_cast<size_t>(tess.GetPointNo()));
      for(int i=0;i<tess.GetPointNo();++i){
	vector<Vector2D> convhull;
	ConvexHull(convhull,&tess,i);
	for(size_t j=0;j<convhull.size();++j){
	  xvert.push_back(convhull[j].x);
	  yvert.push_back(convhull[j].y);
	}
	nvert[static_cast<size_t>(i)]= static_cast<int>(convhull.size());
      }
    }
  };

  class TracerSlice: public Index2Member<double>
  {
  public:

    TracerSlice(const vector<vector<double> >& tracers,
		size_t index):
      tracers_(tracers), index_(index) {}

    size_t getLength(void) const
    {
      return tracers_.size();
    }
		
    double operator()(size_t i) const
    {
      return tracers_[i][index_];
    }

  private:
    const vector<vector<double> >& tracers_;
    const size_t index_;
  };

  vector<double> cold_flows_block(const hdsim& sim)
  {
    vector<double> res(3,0);
    if(sim.GetColdFlowFlag()){
      res[0] = 1;
      sim.GetColdFlowParm(res[1],res[2]);
    }
    else
      res[0] = -1;
    return res;
  }

  vector<double> density_floor_block(const hdsim& sim)
  {
    vector<double> res(3,0);
    if(sim.GetDensityFloorFlag()){
      res[0] = 1;
      sim.GetDensityFloorParm(res[1],res[2]);
    }
    else
      res[0] = -1;
    return res;
  }
}

void write_snapshot_to_hdf5(hdsim const& sim,string const& fname)
{
  ConvexHullData chd(sim.GetTessellation());
  HDF5Shortcut h5sc(fname);
  h5sc("time",vector<double>(1,sim.GetTime()))
    ("x_coordinate",serial_generate
     (MeshGeneratingPointCoordinate(sim.GetTessellation(),&Vector2D::x)))
    ("y_coordinate",serial_generate
     (MeshGeneratingPointCoordinate(sim.GetTessellation(),&Vector2D::y)))
#ifdef RICH_MPI
    ("proc_x_coordinate",serial_generate
     (MeshGeneratingPointCoordinate(sim.GetProcTessellation(),&Vector2D::x)))
    ("proc_y_coordinate", serial_generate
     (MeshGeneratingPointCoordinate(sim.GetProcTessellation(),&Vector2D::y)))
#endif
    ("density",serial_generate
     (CellsPropertyExtractor(sim,ThermalPropertyExtractor(&Primitive::Density))))
    ("pressure",serial_generate
     (CellsPropertyExtractor(sim,ThermalPropertyExtractor(&Primitive::Pressure))))
    ("x_velocity",serial_generate
     (CellsPropertyExtractor(sim,CellVelocityComponentExtractor(&Vector2D::x))))
    ("y_velocity",serial_generate
     (CellsPropertyExtractor(sim,CellVelocityComponentExtractor(&Vector2D::y))))
    ("x position of vertices",chd.xvert)
    ("y position of vertices",chd.yvert)
    ("Number of vertices in cell",chd.nvert)
    ("Cold Flow parameters",cold_flows_block(sim))
    ("Cfl Number",vector<double>(1,sim.GetCfl()))
    ("Cycle number",vector<int>(1,sim.GetCycle()))
    ("Density floor parameters",density_floor_block(sim))
    ("Custom evolution indices",
     list_static_cast<int,size_t>(sim.custom_evolution_indices));

  if(!sim.getTracers().empty()){
    h5sc("Number of tracers",vector<int>(1,static_cast<int>(sim.getTracers()[0].size())));
    for(size_t i=0;i<sim.getTracers()[0].size();++i)
      h5sc("Tracer number "+int2str(static_cast<int>(i)+1),
	   serial_generate(TracerSlice(sim.getTracers(),i)));
  }
  else
    h5sc("Number of tracers",vector<int>(1,0));
}

void read_hdf5_snapshot(ResetDump &dump,string const& fname,EquationOfState
			const* eos)
{
  H5File file(fname, H5F_ACC_RDONLY );
  // Get the mesh points
  vector<double> x=read_double_vector_from_hdf5(file,"x_coordinate");
  vector<double> y=read_double_vector_from_hdf5(file,"y_coordinate");
  dump.snapshot.mesh_points.resize(x.size());
  for(size_t i=0;i<x.size();++i)
    {
      dump.snapshot.mesh_points[i].Set(x[i],y[i]);
    }

  // Get the processor mesh points
#ifdef RICH_MPI
  vector<double> xproc=read_double_vector_from_hdf5(file,"proc_x_coordinate");
  vector<double> yproc=read_double_vector_from_hdf5(file,"proc_y_coordinate");
  dump.procmesh.resize(xproc.size());
  for(int i=0;i<static_cast<int>(xproc.size());++i)
    dump.procmesh[static_cast<size_t>(i)].Set(xproc[static_cast<size_t>(i)],yproc[static_cast<size_t>(i)]);
#endif

  // Get the hydro variables
  vector<double> density=read_double_vector_from_hdf5(file,"density");
  vector<double> pressure=read_double_vector_from_hdf5(file,"pressure");
  vector<double> x_velocity=read_double_vector_from_hdf5(file,"x_velocity");
  vector<double> y_velocity=read_double_vector_from_hdf5(file,"y_velocity");
  dump.snapshot.cells.resize(density.size());
  if(dump.snapshot.cells.size()!=dump.snapshot.mesh_points.size())
    throw UniversalError("Primitive size is not equal to mesh size");
  for(size_t i=0;i<density.size();++i)
    {
      dump.snapshot.cells[i].Density=density[i];
      dump.snapshot.cells[i].Pressure=pressure[i];
      dump.snapshot.cells[i].Velocity=Vector2D(x_velocity[i],y_velocity[i]);
      dump.snapshot.cells[i].Energy=eos->dp2e(dump.snapshot.cells[i].Density,
					      dump.snapshot.cells[i].Pressure);
      dump.snapshot.cells[i].SoundSpeed=eos->dp2c(dump.snapshot.cells[i].Density,
						  dump.snapshot.cells[i].Pressure);
    }
  // read time
  vector<double> time_vector=read_double_vector_from_hdf5(file,"time");
  dump.time=time_vector[0];
  // read tracers if needed
  vector<int> TracerNumber=read_int_vector_from_hdf5(file,"Number of tracers");
  if(TracerNumber[0]>0)
    {
      // resize the tracers
      dump.tracers.resize(x.size());
      for(size_t i=0;i<x.size();++i)
	{
	  dump.tracers[i].resize(static_cast<size_t>(TracerNumber[0]));
	}
      // read the data
      for(int i=0;i<TracerNumber[0];++i)
	{
	  vector<double> tracer=read_double_vector_from_hdf5(file,
							     "Tracer number "+int2str(i+1));
	  if(tracer.size()!=x.size())
	    throw UniversalError("Tracer size not equal to mesh size");
	  for(size_t j=0;j<x.size();++j)
	    dump.tracers[j][static_cast<size_t>(i)]=tracer[j];
	}
    }
  // read the coldflows parameters
  vector<double> coldflows=read_double_vector_from_hdf5(file,
							"Cold Flow parameters");
  if(coldflows[0]>0)
    dump.coldflows=true;
  else
    dump.coldflows=false;
  dump.a=coldflows[1];
  dump.b=coldflows[2];
  // read the cfl
  vector<double> cfl=read_double_vector_from_hdf5(file,"Cfl Number");
  dump.cfl=cfl[0];
  // read the cycle number
  vector<int> cycle_number=read_int_vector_from_hdf5(file,"Cycle number");
  dump.cycle=cycle_number[0];
  // read the density floor parameters
  vector<double> densityfloor=read_double_vector_from_hdf5(file,
							   "Density floor parameters");
  if(densityfloor[0]>0)
    dump.densityfloor=true;
  else
    dump.densityfloor=false;
  dump.densitymin=densityfloor[1];
  dump.pressuremin=densityfloor[2];
  // read the custom evolution indeces
  dump.cevolve = read_sizet_vector_from_hdf5(file, "Custom evolution indices");
}

void ConvertHDF5toBinary(string const& input, string const& output)
{
	H5File file(input, H5F_ACC_RDONLY);
	fstream myFile(output.c_str(), ios::out | ios::binary);
	// Get the mesh points
	vector<double> x = read_double_vector_from_hdf5(file, "x_coordinate");
	vector<double> y = read_double_vector_from_hdf5(file, "y_coordinate");
	int temp = static_cast<int>(x.size());
	myFile.write(reinterpret_cast<char*>(&temp), sizeof (int));
	for (int i = 0; i<temp; ++i)
	{
	  myFile.write(reinterpret_cast<char*>(&x[static_cast<size_t>(i)]), sizeof(double));
	  myFile.write(reinterpret_cast<char*>(&y[static_cast<size_t>(i)]), sizeof(double));
	}
	// Get the processor mesh points
#ifdef RICH_MPI
	vector<double> xproc = read_double_vector_from_hdf5(file, "proc_x_coordinate");
	vector<double> yproc = read_double_vector_from_hdf5(file, "proc_y_coordinate");
	int temp2 = static_cast<int>(xproc.size());
	myFile.write(reinterpret_cast<char*>(&temp2), sizeof (int));
	for (int i = 0; i<temp2; ++i)
	{
	  myFile.write(reinterpret_cast<char*>(&x[static_cast<size_t>(i)]), sizeof(double));
	  myFile.write(reinterpret_cast<char*>(&y[static_cast<size_t>(i)]), sizeof(double));
	}
#endif

	// Get the hydro variables
	vector<double> density = read_double_vector_from_hdf5(file, "density");
	vector<double> pressure = read_double_vector_from_hdf5(file, "pressure");
	vector<double> x_velocity = read_double_vector_from_hdf5(file, "x_velocity");
	vector<double> y_velocity = read_double_vector_from_hdf5(file, "y_velocity");
	for (int i = 0; i<temp; ++i)
	{
	  myFile.write(reinterpret_cast<char*>(&pressure[static_cast<size_t>(i)]), sizeof(double));
	  myFile.write(reinterpret_cast<char*>(&density[static_cast<size_t>(i)]), sizeof(double));
	  myFile.write(reinterpret_cast<char*>(&x_velocity[static_cast<size_t>(i)]), sizeof(double));
	  myFile.write(reinterpret_cast<char*>(&y_velocity[static_cast<size_t>(i)]), sizeof(double));
	}
	// read time
	vector<double> time_vector = read_double_vector_from_hdf5(file, "time");
	myFile.write(reinterpret_cast<char*>(&time_vector[0]), sizeof(double));
	// read the coldflows parameters
	vector<double> coldflows = read_double_vector_from_hdf5(file,"Cold Flow parameters");
	char cold = static_cast<char>(coldflows[0]);
	myFile.write(reinterpret_cast<char*>(&cold), sizeof(char));
	// read the cfl
	vector<double> cfl = read_double_vector_from_hdf5(file, "Cfl number");
	myFile.write(reinterpret_cast<char*>(&cfl[0]), sizeof(double));
	myFile.write(reinterpret_cast<char*>(&coldflows[1]), sizeof(double));
	myFile.write(reinterpret_cast<char*>(&coldflows[2]), sizeof(double));
	// read the cycle number
	vector<int> cycle_number = read_int_vector_from_hdf5(file, "Cycle number");
	myFile.write(reinterpret_cast<char*>(&cycle_number[0]), sizeof(int));
	// read tracers if needed
	vector<int> TracerNumber = read_int_vector_from_hdf5(file, "Number of tracers");
	myFile.write(reinterpret_cast<char*>(&TracerNumber[0]), sizeof(int));
	// read the density floor parameters
	vector<double> densityfloor = read_double_vector_from_hdf5(file,"Density floor parameters");
	cold = (densityfloor[0]>0) ? '1' : '0';
	myFile.write(reinterpret_cast<char*>(&cold), sizeof(char));
	myFile.write(reinterpret_cast<char*>(&densityfloor[1]), sizeof(double));
	myFile.write(reinterpret_cast<char*>(&densityfloor[2]), sizeof(double));
	// read the custom evolution indeces
	vector<size_t> cevolve = read_sizet_vector_from_hdf5(file, "Custom evolution indeces");
	for (int i = 0; i<temp; ++i)
	{
	  myFile.write(reinterpret_cast<char*>(&cevolve[static_cast<size_t>(i)]), sizeof(unsigned int));
	}
	if (TracerNumber[0]>0)
	{
		// read the data
		for (int i = 0; i<TracerNumber[0]; ++i)
		{
			vector<double> tracer = read_double_vector_from_hdf5(file, "Tracer number " + int2str(i + 1));
			if (tracer.size() != x.size())
				throw UniversalError("Tracer size not equal to mesh size");
			for (size_t j = 0; j<tracer.size(); ++j)
			  myFile.write(reinterpret_cast<char*>(&tracer[j]), sizeof(double));
		}
	}
	myFile.close();
}
