#include "hdf5_diagnostics.hpp"
#include "../../misc/hdf5_utils.hpp"
#include "../../misc/lazy_list.hpp"

using namespace H5;

DiagnosticAppendix::~DiagnosticAppendix(void) {}

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
  class MeshGeneratingPointCoordinate: public LazyList<double>
  {
  public:

    MeshGeneratingPointCoordinate(const Tessellation& tess,
				  double Vector2D::* component):
      tess_(tess), component_(component) {}

    size_t size(void) const
    {
      return static_cast<size_t>(tess_.GetPointNo());
    }

    double operator[](size_t i) const
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

    virtual double operator()(const ComputationalCell& p) const = 0;

    virtual ~SingleCellPropertyExtractor(void) {}
  };

  class ThermalPropertyExtractor: public SingleCellPropertyExtractor
  {
  public:

    explicit ThermalPropertyExtractor(double ComputationalCell::* var):
      var_(var) {}

    double operator()(const ComputationalCell& p) const
    {
      return p.*var_;
    }

  private:
    double ComputationalCell::* var_;
  };

  class CellVelocityComponentExtractor: public SingleCellPropertyExtractor
  {
  public:

    explicit CellVelocityComponentExtractor(double Vector2D::* component):
      component_(component) {}

    double operator()(const ComputationalCell& p) const
    {
      return p.velocity.*component_;
    }

  private:
    double Vector2D::* component_;
  };

  class CellsPropertyExtractor: public LazyList<double>
  {
  public:

    CellsPropertyExtractor(const hdsim& sim,
			   const SingleCellPropertyExtractor& scpe):
      sim_(sim), scpe_(scpe) {}

    size_t size(void) const
    {
      return static_cast<size_t>(sim_.getTessellation().GetPointNo());
    }

    double operator[](size_t i) const
    {
      return scpe_(sim_.getAllCells()[i]);
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

    explicit ConvexHullData(const Tessellation& tess):
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

  class StickerSlice: public LazyList<double>
  {
  public:

    StickerSlice(const hdsim& sim,
		 const string& name):
      sim_(sim), name_(name) {}

    size_t size(void) const
    {
      return sim_.getAllCells().size();
    }

    double operator[](size_t i) const
    {
      return static_cast<double>
	(safe_retrieve
	 (sim_.getAllCells()[i].stickers,name_));
    }

  private:
    const hdsim& sim_;
    const string& name_;
  };

  class TracerSlice: public LazyList<double>
  {
  public:

    TracerSlice(const hdsim& sim,
		const string& name):
      sim_(sim), name_(name) {}

    size_t size(void) const
    {
      return sim_.getAllCells().size();
    }
		
    double operator[](size_t i) const
    {
      return safe_retrieve
	(sim_.getAllCells()[i].tracers,
	 name_);
    }

  private:
    const hdsim& sim_;
    const string& name_;
  };
}

void write_snapshot_to_hdf5(hdsim const& sim,string const& fname,
			    const vector<DiagnosticAppendix*>& appendices)
{
  ConvexHullData chd(sim.getTessellation());
  HDF5Shortcut h5sc(fname);
  h5sc("time",vector<double>(1,sim.getTime()))
    ("x_coordinate",serial_generate
     (MeshGeneratingPointCoordinate(sim.getTessellation(),&Vector2D::x)))
    ("y_coordinate",serial_generate
     (MeshGeneratingPointCoordinate(sim.getTessellation(),&Vector2D::y)))
#ifdef RICH_MPI
    ("proc_x_coordinate",serial_generate
     (MeshGeneratingPointCoordinate(sim.GetProcTessellation(),&Vector2D::x)))
    ("proc_y_coordinate", serial_generate
     (MeshGeneratingPointCoordinate(sim.GetProcTessellation(),&Vector2D::y)))
#endif
    ("density",serial_generate
     (CellsPropertyExtractor
      (sim,ThermalPropertyExtractor(&ComputationalCell::density))))
    ("pressure",serial_generate
     (CellsPropertyExtractor
      (sim,ThermalPropertyExtractor(&ComputationalCell::pressure))))
    ("x_velocity",serial_generate
     (CellsPropertyExtractor(sim,CellVelocityComponentExtractor(&Vector2D::x))))
    ("y_velocity",serial_generate
     (CellsPropertyExtractor(sim,CellVelocityComponentExtractor(&Vector2D::y))))
    ("x position of vertices",chd.xvert)
    ("y position of vertices",chd.yvert)
    ("Number of vertices in cell",chd.nvert)
    ("Cycle number",vector<int>(1,sim.getCycle()));

  h5sc("Number of tracers", vector<int>
       (1,
	static_cast<int>(sim.getAllCells().front().tracers.size())));

  for(boost::container::flat_map<std::string,double>::const_iterator it=
	sim.getAllCells().front().tracers.begin();
      it!=sim.getAllCells().front().tracers.end(); ++it)
    h5sc(it->first,serial_generate(TracerSlice(sim,it->first)));

  for(boost::container::flat_map<std::string,bool>::const_iterator it=
	sim.getAllCells().front().stickers.begin();
      it!=sim.getAllCells().front().stickers.end(); ++it)
    h5sc(it->first,serial_generate(StickerSlice(sim,it->first)));

  for(size_t i=0;i<appendices.size();++i)
    h5sc(appendices[i]->getName(),(*(appendices[i]))(sim));
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
