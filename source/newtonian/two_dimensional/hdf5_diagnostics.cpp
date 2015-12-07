#include "hdf5_diagnostics.hpp"
#include "../../misc/hdf5_utils.hpp"
#include "../../misc/lazy_list.hpp"
#ifdef RICH_MPI
#include "../../tessellation/VoronoiMesh.hpp"
#include <boost/mpi/communicator.hpp>
#endif

using namespace H5;

Snapshot::Snapshot(void):
  mesh_points(),
  cells(),
  time(),
  cycle(){}

Snapshot::Snapshot(const Snapshot& source):
  mesh_points(source.mesh_points),
  cells(source.cells),
  time(source.time),
  cycle(source.cycle){}

DiagnosticAppendix::~DiagnosticAppendix(void) {}

namespace {

  template<class T> vector<T> read_vector_from_hdf5
  (const CommonFG& file,
   const string& caption,
   const DataType& datatype)
  {
    DataSet dataset = file.openDataSet(caption);
    DataSpace filespace = dataset.getSpace();
    hsize_t dims_out[2];
    filespace.getSimpleExtentDims(dims_out,NULL);
    const size_t NX = static_cast<size_t>(dims_out[0]);
    vector<T> result(NX);
    dataset.read(&result[0],datatype);
    return result;
  }

  vector<double> read_double_vector_from_hdf5
  (CommonFG& file,string const& caption)
  {
    return read_vector_from_hdf5<double>
      (file,
       caption,
       PredType::NATIVE_DOUBLE);
  }

  vector<int> read_int_vector_from_hdf5
  (const CommonFG& file,
   const string& caption)
  {
    return read_vector_from_hdf5<int>
      (file,
       caption,
       PredType::NATIVE_INT);
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
	ConvexHull(convhull,tess,i);
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
		return static_cast<size_t>(sim_.getTessellation().GetPointNo());
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
      return static_cast<size_t>(sim_.getTessellation().GetPointNo());
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
  H5File file(H5std_string(fname), H5F_ACC_TRUNC);
  Group geometry = file.createGroup("/geometry");
  Group gappendices = file.createGroup("/appendices");
  Group hydrodynamic = file.createGroup("/hydrodynamic");
  Group tracers = file.createGroup("/tracers");
  Group stickers = file.createGroup("/stickers");
#ifdef RICH_MPI
  Group mpi = file.createGroup("/mpi");
#endif

  // General
  write_std_vector_to_hdf5
    (file,
     vector<double>(1,sim.getTime()),
     "time");
  write_std_vector_to_hdf5
    (file,
     vector<int>(1,sim.getCycle()),
     "cycle");

  // Geometry  
  write_std_vector_to_hdf5
    (geometry,
     serial_generate
     (MeshGeneratingPointCoordinate
      (sim.getTessellation(),&Vector2D::x)),
     "x_coordinate");
  write_std_vector_to_hdf5
    (geometry,
     serial_generate
     (MeshGeneratingPointCoordinate
      (sim.getTessellation(),&Vector2D::y)),
     "y_coordinate");
    write_std_vector_to_hdf5
    (geometry,
     chd.xvert,
     "x_vertices");
  write_std_vector_to_hdf5
    (geometry,
     chd.yvert,
     "y_vertices");
  write_std_vector_to_hdf5
    (geometry,
     chd.nvert,
     "n_vertices");
  //MPI
#ifdef RICH_MPI
  write_std_vector_to_hdf5
	  (mpi,
		  serial_generate
		  (MeshGeneratingPointCoordinate
			  (sim.GetProcTessellation(), &Vector2D::x)),
		  "x_coordinate");
  write_std_vector_to_hdf5
	  (mpi,
		  serial_generate
		  (MeshGeneratingPointCoordinate
			  (sim.GetProcTessellation(), &Vector2D::y)),
		  "y_coordinate");
#endif

  // Hydrodynamic
  write_std_vector_to_hdf5
    (hydrodynamic,
     serial_generate
     (CellsPropertyExtractor
      (sim,ThermalPropertyExtractor(&ComputationalCell::density))),
     "density");
  write_std_vector_to_hdf5
    (hydrodynamic,
     serial_generate
     (CellsPropertyExtractor
      (sim,ThermalPropertyExtractor(&ComputationalCell::pressure))),
     "pressure");
  write_std_vector_to_hdf5
    (hydrodynamic,
     serial_generate
     (CellsPropertyExtractor
      (sim,CellVelocityComponentExtractor(&Vector2D::x))),
     "x_velocity");
  write_std_vector_to_hdf5
    (hydrodynamic,
     serial_generate
     (CellsPropertyExtractor
      (sim,CellVelocityComponentExtractor(&Vector2D::y))),
     "y_velocity");

  // Tracers
  for(boost::container::flat_map<std::string,double>::const_iterator it=
	sim.getAllCells().front().tracers.begin();
      it!=sim.getAllCells().front().tracers.end(); ++it)
    write_std_vector_to_hdf5
      (tracers,
       serial_generate(TracerSlice(sim,it->first)),
       it->first);

  // Stickers
  for(boost::container::flat_map<std::string,bool>::const_iterator it=
	sim.getAllCells().front().stickers.begin();
      it!=sim.getAllCells().front().stickers.end(); ++it)
    write_std_vector_to_hdf5
      (stickers,
       serial_generate(StickerSlice(sim,it->first)),
       it->first);

  // Appendices
  for(size_t i=0;i<appendices.size();++i)
    write_std_vector_to_hdf5
      (gappendices,
       (*(appendices.at(i)))(sim),
       appendices.at(i)->getName());
}

Snapshot read_hdf5_snapshot
(const string& fname)
{
  Snapshot res;
  H5File file(fname, H5F_ACC_RDONLY );
  Group g_geometry = file.openGroup("geometry");
  Group g_hydrodynamic = file.openGroup("hydrodynamic");
  Group g_tracers = file.openGroup("tracers");
  Group g_stickers = file.openGroup("stickers");
#ifdef RICH_MPI
  Group mpi = file.createGroup("/mpi");
#endif

  
  // Mesh points
  {
    const vector<double> x=
      read_double_vector_from_hdf5(g_geometry,"x_coordinate");
    const vector<double> y=
      read_double_vector_from_hdf5(g_geometry,"y_coordinate");
    res.mesh_points.resize(x.size());
    for(size_t i=0;i<x.size();++i)
      res.mesh_points.at(i) = Vector2D(x.at(i), y.at(i));
  }

#ifdef RICH_MPI
  // MPI
  {
	  const vector<double> x =
		  read_double_vector_from_hdf5(mpi, "x_coordinate");
	  const vector<double> y =
		  read_double_vector_from_hdf5(mpi, "y_coordinate");
	  res.proc_points.resize(x.size());
	  for (size_t i = 0; i<x.size(); ++i)
		  res.proc_points.at(i) = Vector2D(x.at(i), y.at(i));
  }
#endif

  // Hydrodynamic
  {
    const vector<double> density =
      read_double_vector_from_hdf5(g_hydrodynamic,"density");
    const vector<double> pressure =
      read_double_vector_from_hdf5(g_hydrodynamic,"pressure");
    const vector<double> x_velocity =
      read_double_vector_from_hdf5(g_hydrodynamic,"x_velocity");
    const vector<double> y_velocity =
      read_double_vector_from_hdf5(g_hydrodynamic,"y_velocity");
    boost::container::flat_map<string,vector<double> > tracers;
    for(hsize_t n=0;n<g_tracers.getNumObjs();++n){
      const H5std_string name = g_tracers.getObjnameByIdx(n);
      tracers[name] = read_double_vector_from_hdf5
	(g_tracers,name);
    }
    boost::container::flat_map<string,vector<int> > stickers;
    for(hsize_t n=0;n<g_stickers.getNumObjs();++n){
      const H5std_string name = g_stickers.getObjnameByIdx(n);
      stickers[name] =
	read_int_vector_from_hdf5(g_stickers,name);
    }
    res.cells.resize(density.size());
    for(size_t i=0;i<res.cells.size();++i){
      res.cells.at(i).density = density.at(i);
      res.cells.at(i).pressure = pressure.at(i);
      res.cells.at(i).velocity.x = x_velocity.at(i);
      res.cells.at(i).velocity.y = y_velocity.at(i);
      for(boost::container::flat_map<string,vector<double> >::const_iterator it = tracers.begin();
	  it!=tracers.end();
	  ++it)
	res.cells.at(i).tracers[it->first] = (it->second).at(i);
      for(boost::container::flat_map<string,vector<int> >::const_iterator it = stickers.begin();
	  it!=stickers.end();
	  ++it)
	res.cells.at(i).stickers[it->first] = ((it->second).at(i) == 1);
    }
  }

  // Misc
  {
    const vector<double> time =
      read_double_vector_from_hdf5(file,"time");
    res.time = time.at(0);
	const vector<int> cycle =
		read_int_vector_from_hdf5(file, "cycle");
	res.cycle = cycle.at(0);
  }
    
  return res;
}

void WriteDelaunay(Delaunay const& tri, string const& filename)
{
	vector<Vector2D> const& cor=tri.getCor();
	vector<double> x_cor, y_cor;
	vector<int> facets;
	int nfacets=tri.get_num_facet();

	H5File file(H5std_string(filename), H5F_ACC_TRUNC);

	for (size_t i = 0; i < cor.size(); ++i)
	{
		x_cor.push_back(cor[i].x);
		y_cor.push_back(cor[i].y);
	}

	for (int i = 0; i < nfacets; ++i)
	{
		facets.push_back(tri.get_facet(i).vertices.first);
		facets.push_back(tri.get_facet(i).vertices.second);
		facets.push_back(tri.get_facet(i).vertices.third);
	}

	write_std_vector_to_hdf5(file, x_cor, "x_coordinate");
	write_std_vector_to_hdf5(file, y_cor, "y_coordinate");
	write_std_vector_to_hdf5(file, vector<int>(1, tri.GetOriginalLength()), "point number");
	write_std_vector_to_hdf5(file, facets, "triangles");
}

#ifdef RICH_MPI

namespace
{

}

Snapshot ReDistributeData(string const& filename, Tessellation const& proctess,size_t snapshot_number)
{
	const boost::mpi::communicator world;
	double read_num = 
	  static_cast<double>(snapshot_number)*1.0/
	  static_cast<double>(world.size());
	// Read the data
	int start = static_cast<int>
	  (floor
	   (static_cast<double>(world.rank())*read_num+0.1));
	int stop = static_cast<int>
	  (floor((1+world.rank())*read_num-1.1));
	Snapshot res,snap;
	for (int i = start; i < stop; ++i)
	{
		Snapshot temp = read_hdf5_snapshot(filename + int2str(i) + ".h5");
		snap.cells.insert(snap.cells.end(), temp.cells.begin(), temp.cells.end());
		snap.mesh_points.insert(snap.mesh_points.end(), temp.mesh_points.begin(), temp.mesh_points.end());
		if (i == start)
		{
			snap.time = temp.time;
			snap.cycle = temp.cycle;
		}
	}
	vector<vector<Vector2D> > chull(static_cast<size_t>(world.size()));
	vector<vector<ComputationalCell> > cell_recv(world.size());
	vector<vector<Vector2D> > mesh_recv(world.size());
	vector<vector<size_t> > indeces(world.size());
	for (size_t i = 0; i < chull.size(); ++i)
		ConvexHull(chull[i], proctess, static_cast<int>(i));
	for (size_t i = 0; i < snap.mesh_points.size(); ++i)
	{
		bool added = false;
		for (size_t j = 0; j < chull.size(); ++j)
		{
			if (PointInCell(chull[j], snap.mesh_points[i]))
			{
				indeces[j].push_back(i);
				added = true;
				break;
			}
		}
		if (!added)
			throw UniversalError("Didn't find point in ReDistributeData");
	}
	// Send/Recv data
	vector < boost::mpi::request> req;
	for (size_t i = 0; i < chull.size(); ++i)
	{
		if (i == static_cast<size_t>(world.rank()))
		{ 
			res.cells = VectorValues(snap.cells, indeces[i]);
			res.mesh_points = VectorValues(snap.mesh_points, indeces[i]);
			continue;
		}	
		req.push_back(world.isend(static_cast<int>(i), 0, VectorValues(snap.cells,indeces[i])));
		req.push_back(world.isend(static_cast<int>(i), 1, VectorValues(snap.mesh_points, indeces[i])));
		req.push_back(world.irecv(static_cast<int>(i), 0,cell_recv[i]));
		req.push_back(world.irecv(static_cast<int>(i), 1, mesh_recv[i]));
	}
	boost::mpi::wait_all(req.begin(), req.end());

	for (size_t i = 0; i < chull.size(); ++i)
	{
		if (i == static_cast<size_t>(world.rank()))
			continue;
		res.cells.insert(res.cells.end(), cell_recv[i].begin(), cell_recv[i].end());
		res.mesh_points.insert(res.mesh_points.end(), mesh_recv[i].begin(), mesh_recv[i].end());
	}

	res.time = snap.time;
	res.cycle = snap.cycle;
	return res;
}
#endif

void WriteTess(Tessellation const& tess, string const& fname)
{
	ConvexHullData chd(tess);
	H5File file(H5std_string(fname), H5F_ACC_TRUNC);
	Group geometry = file.createGroup("/geometry");
	write_std_vector_to_hdf5(geometry,serial_generate(MeshGeneratingPointCoordinate(tess, &Vector2D::x)),
		"x_coordinate");
	write_std_vector_to_hdf5(geometry,serial_generate(MeshGeneratingPointCoordinate(tess, &Vector2D::y)),
		"y_coordinate");
	write_std_vector_to_hdf5(geometry,chd.xvert,"x_vertices");
	write_std_vector_to_hdf5(geometry,chd.yvert,"y_vertices");
	write_std_vector_to_hdf5(geometry,chd.nvert,"n_vertices");
}
