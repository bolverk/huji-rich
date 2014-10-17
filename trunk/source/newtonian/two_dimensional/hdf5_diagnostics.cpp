#include "hdf5_diagnostics.hpp"

using namespace H5;

namespace {
	void write_std_vector_to_hdf5
  (H5File& file,
   vector<size_t> const& num_list,
   string const& caption)
  {
    hsize_t dimsf[1];
    dimsf[0] = static_cast<hsize_t>(num_list.size());
    DataSpace dataspace(1, dimsf);

    FloatType datatype(PredType::NATIVE_UINT);
    datatype.setOrder(H5T_ORDER_LE);

    // Modify dataset creation property to enable chunking
    DSetCreatPropList  plist;
    if(dimsf[0]>100000)
      dimsf[0]=100000;
    plist.setChunk(1,dimsf);
    plist.setDeflate(6);

    DataSet dataset = file.createDataSet(H5std_string(caption),
					 datatype,
					 dataspace,plist);

    const vector<unsigned> buf = list_static_cast<unsigned,size_t>(num_list);
    dataset.write(&buf[0],PredType::NATIVE_UINT);
  }

  void write_std_vector_to_hdf5
  (H5File& file,
   vector<double> const& num_list,
   string const& caption)
  {
    hsize_t dimsf[1];
    dimsf[0] = (hsize_t)num_list.size();
    DataSpace dataspace(1, dimsf);

    FloatType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);

    // Modify dataset creation property to enable chunking
    DSetCreatPropList  plist;
    if(dimsf[0]>100000)
      dimsf[0]=100000;
    plist.setChunk(1,dimsf);
    plist.setDeflate(6);

    DataSet dataset = file.createDataSet(H5std_string(caption),
					 datatype,
					 dataspace,plist);

    dataset.write(&num_list[0],PredType::NATIVE_DOUBLE);
  }

  vector<double> read_double_vector_from_hdf5(H5File& file,string const& data_name)
  {
    DataSet dataset = file.openDataSet(data_name);
    DataSpace filespace = dataset.getSpace();
    hsize_t dims_out[2];
    filespace.getSimpleExtentDims(dims_out,NULL);
    int NX = (int)dims_out[0];
    vector<double> result((size_t)NX);
    dataset.read(&result[0],PredType::NATIVE_DOUBLE);
    return result;
  }

  vector<size_t> read_sizet_vector_from_hdf5(H5File& file,string const& data_name)
  {
    DataSet dataset = file.openDataSet(data_name);
    DataSpace filespace = dataset.getSpace();
    hsize_t dims_out[2];
    filespace.getSimpleExtentDims(dims_out,NULL);
    int NX = (int)dims_out[0];
    vector<unsigned> result((size_t)NX);
    dataset.read(&result[0],PredType::NATIVE_UINT);
    return list_static_cast<size_t,unsigned>(result);
  }

  vector<int> read_int_vector_from_hdf5(H5File& file,string const& data_name)
  {
    DataSet dataset = file.openDataSet(data_name);
    DataSpace filespace = dataset.getSpace();
    hsize_t dims_out[2];
    filespace.getSimpleExtentDims(dims_out,NULL);
    int NX=(int)dims_out[0];
    vector<int> result((size_t)NX);
    dataset.read(&result[0],PredType::NATIVE_INT);
    return result;
  }

  void write_std_vector_to_hdf5
  (H5File& file,
   vector<int> const& num_list,
   string const& caption)
  {
    hsize_t dimsf[1];
    dimsf[0] = (hsize_t)num_list.size();
    DataSpace dataspace(1, dimsf);

    IntType datatype(PredType::NATIVE_INT);
    datatype.setOrder(H5T_ORDER_LE);

    // Modify dataset creation property to enable chunking
    DSetCreatPropList  plist;
    if(dimsf[0]>100000)
      dimsf[0]=100000;
    plist.setChunk(1,dimsf);
    plist.setDeflate(6);

    DataSet dataset = file.createDataSet(H5std_string(caption),
					 datatype,
					 dataspace,plist);

    dataset.write(&num_list[0],PredType::NATIVE_INT);
  }
}

void write_snapshot_to_hdf5(hdsim const& sim,string const& fname)
{
  // Create file
  H5File file(H5std_string(fname),
	      H5F_ACC_TRUNC);

  // Write time
  {
    vector<double> time_vector(1,0);
    time_vector[0] = sim.GetTime();
    write_std_vector_to_hdf5(file, time_vector, "time");
  }

  // Write grid
  {
    vector<double> x_coordinate(sim.GetCellNo(),0);
    vector<double> y_coordinate(sim.GetCellNo(),0);
    for(int i=0;i<sim.GetCellNo();++i){
      x_coordinate[(size_t)i] = sim.GetMeshPoint(i).x;
      y_coordinate[(size_t)i] = sim.GetMeshPoint(i).y;
    }
    write_std_vector_to_hdf5(file, x_coordinate, "x_coordinate");
    write_std_vector_to_hdf5(file, y_coordinate, "y_coordinate");
  }

  // write processor grid
#ifdef RICH_MPI
  const int nproc=sim.GetProcTessellation().GetPointNo();
  vector<double> xproc(nproc),yproc(nproc);
  for(int i=0;i<nproc;++i)
  {
	  xproc[i]=sim.GetProcTessellation().GetMeshPoint(i).x;
	  yproc[i]=sim.GetProcTessellation().GetMeshPoint(i).y;
  }
  write_std_vector_to_hdf5(file,xproc, "proc_x_coordinate");
  write_std_vector_to_hdf5(file,yproc, "proc_y_coordinate");
#endif

  // Write hydrodynamic variables
  {
    vector<double> density_list((size_t)sim.GetCellNo());
    vector<double> pressure_list((size_t)sim.GetCellNo());
    vector<double> x_velocity_list((size_t)sim.GetCellNo());
    vector<double> y_velocity_list((size_t)sim.GetCellNo());
    for(int i=0;i<sim.GetCellNo();++i){
      density_list[(size_t)i] = sim.GetCell(i).Density;
      pressure_list[(size_t)i] = sim.GetCell(i).Pressure;
      x_velocity_list[(size_t)i] = sim.GetCell(i).Velocity.x;
      y_velocity_list[(size_t)i] = sim.GetCell(i).Velocity.y;
    }
    write_std_vector_to_hdf5(file,density_list,"density");
    write_std_vector_to_hdf5(file,pressure_list,"pressure");
    write_std_vector_to_hdf5(file,x_velocity_list,"x_velocity");
    write_std_vector_to_hdf5(file,y_velocity_list,"y_velocity");
  }
  // write cell vertices
  // Do the convex hull for each point
  vector<Vector2D> convhull;
  vector<double> xvert,yvert;
  vector<int> nvert((size_t)sim.GetCellNo());
  xvert.reserve((size_t)(7*sim.GetCellNo()));
  yvert.reserve((size_t)(7*sim.GetCellNo()));
  for(int i=0;i<sim.GetCellNo();++i)
    {
      ConvexHull(convhull,&sim.GetTessellation(),i);
      for(int j=0;j<(int)convhull.size();++j)
	{
	  xvert.push_back(convhull[(size_t)j].x);
	  yvert.push_back(convhull[(size_t)j].y);
	}
      nvert[(size_t)i]=(int)convhull.size();
    }
  write_std_vector_to_hdf5(file,xvert,"x position of vertices");
  write_std_vector_to_hdf5(file,yvert,"y position of vertices");
  write_std_vector_to_hdf5(file,nvert,"Number of vertices in cell");

  // write tracers if needed
  vector<vector<double> > tracers=sim.getTracers();
  vector<int> number_of_tracers;
  number_of_tracers.push_back(0);
  if(!tracers.empty())
    number_of_tracers[0]=(int)tracers[0].size();
  write_std_vector_to_hdf5(file,number_of_tracers,"Number of tracers");
  if(!tracers.empty())
    {
      for(int i=0;i<(int)tracers[0].size();++i)
	{
	  vector<double> tracer_temp(tracers.size(),0);
	  for(int j=0;j<(int)tracers.size();++j)
	    tracer_temp[(size_t)j]=tracers[(size_t)j][(size_t)i];
	  write_std_vector_to_hdf5(file,tracer_temp,"Tracer number "+int2str(i+1));
	}
    }
  // write the coldflows parameters
  vector<double> coldflows(3,0);
  if(sim.GetColdFlowFlag())
    {
      coldflows[0]=1;
      double a,b;
      sim.GetColdFlowParm(a,b);
      coldflows[1]=a;
      coldflows[2]=b;
    }
  else
    coldflows[0]=-1;
  write_std_vector_to_hdf5(file,coldflows,"Cold Flow parameters");
  // write the cfl
  vector<double> cfl(1,0);
  cfl[0]=sim.GetCfl();
  write_std_vector_to_hdf5(file,cfl,"Cfl number");
  // write the cycle number
  vector<int> cycle_number(1,0);
  cycle_number[0]=sim.GetCycle();
  write_std_vector_to_hdf5(file,cycle_number,"Cycle number");
  // write the density floor parameters
  vector<double> densityfloor(3,0);
  if(sim.GetDensityFloorFlag())
    {
      densityfloor[0]=1;
      double d_min,p_min;
      sim.GetDensityFloorParm(d_min,p_min);
      densityfloor[1]=d_min;
      densityfloor[2]=p_min;
    }
  else
    densityfloor[0]=-1;
  write_std_vector_to_hdf5(file,densityfloor,"Density floor parameters");
  // write the custom evolution indeces
  write_std_vector_to_hdf5(file,sim.custom_evolution_indices,"Custom evolution indeces");
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
  for(int i=0;i<(int)xproc.size();++i)
      dump.procmesh[i].Set(xproc[i],yproc[i]);
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
	  dump.tracers[i].resize((size_t)TracerNumber[0]);
	}
      // read the data
      for(int i=0;i<TracerNumber[0];++i)
	{
	  vector<double> tracer=read_double_vector_from_hdf5(file,
							     "Tracer number "+int2str(i+1));
	  if(tracer.size()!=x.size())
	    throw UniversalError("Tracer size not equal to mesh size");
	  for(size_t j=0;j<x.size();++j)
	    dump.tracers[j][(size_t)i]=tracer[j];
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
  vector<double> cfl=read_double_vector_from_hdf5(file,"Cfl number");
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
  dump.cevolve= read_sizet_vector_from_hdf5(file,"Custom evolution indeces");
}
