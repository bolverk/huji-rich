#include "hdf5_diagnostics1d.hpp"

using H5::PredType;
using H5::DataSet;
using H5::FloatType;
using H5::DataSpace;
using H5::H5File;
using std::vector;

namespace {
  void write_std_vector_to_hdf5
  (H5File& file,
   vector<double> const& num_list,
   string const& caption)
  {
    hsize_t dimsf[1];
    dimsf[0] = num_list.size();
    DataSpace dataspace(1, dimsf);
    
    FloatType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);

    DataSet dataset = file.createDataSet(H5std_string(caption),
					 datatype,
					 dataspace);

    /* To do: Replace raw pointers with smart pointers
     */
    double *data = new double[static_cast<int>(num_list.size())];
    for(size_t i=0;i<num_list.size();++i)
      data[i] = num_list[i];
    dataset.write(data, PredType::NATIVE_DOUBLE);
    delete[] data;
  }
}

void diagnostics1d::write_snapshot_to_hdf5
(hdsim1D const& sim,
 string const& fname)
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
    const vector<double>& grid = sim.getState().getVertices();
    vector<double> grid_vector(grid.size()-1);
    for(size_t i=0;i<grid_vector.size();++i)
      grid_vector[i] = 0.5*(grid.at(i)+grid.at(i+1));
    write_std_vector_to_hdf5(file, grid_vector, "grid");
  }

  // Write Hydrodynamic variables
  {
    const vector<ComputationalCell>& cells = sim.getState().getCells();
    const size_t n = cells.size();
    vector<double> density_vector(n);
    vector<double> pressure_vector(n);
    vector<double> x_velocity_vector(n);
    vector<double> y_velocity_vector(n);
    //    const vector<Primitive>& cells = sim.getCells();
    for(size_t i=0;i<n;++i){
      density_vector[size_t(i)] = cells.at(i).density;
      pressure_vector[size_t(i)] = cells.at(i).pressure;
      x_velocity_vector[size_t(i)] = cells.at(i).velocity.x;
      y_velocity_vector[size_t(i)] = cells.at(i).velocity.y;
    }
    write_std_vector_to_hdf5(file, density_vector, "density");
    write_std_vector_to_hdf5(file, pressure_vector, "pressure");
    write_std_vector_to_hdf5(file, x_velocity_vector, "x_velocity");
    write_std_vector_to_hdf5(file, y_velocity_vector, "y_velocity");
  }

  // Write tracers
  /*
  {
    const size_t n = sim.getCells().size();
    const size_t m = sim.getS
  }
  */
}
