#include "hdf5_diagnostics1d.hpp"

namespace {
  void write_std_vector_to_hdf5
  (H5File& file,
   vector<double> const& num_list,
   string const& caption)
  {
    hsize_t dimsf[1];
    dimsf[0] = (int)num_list.size();
    DataSpace dataspace(1, dimsf);
    
    FloatType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);

    DataSet dataset = file.createDataSet(H5std_string(caption),
					 datatype,
					 dataspace);

    /* To do: Replace raw pointers with smart pointers
     */
    double *data = new double[(int)num_list.size()];
    for(int i=0;i<(int)num_list.size();++i)
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
    vector<double> grid_vector(sim.GetCellNo());
    for(int i=0;i<sim.GetCellNo();++i)
      grid_vector[i] = sim.GetCellCenter(i);
    write_std_vector_to_hdf5(file, grid_vector, "grid");
  }

  // Write Hydrodynamic variables
  {
    vector<double> density_vector(sim.GetCellNo());
    vector<double> pressure_vector(sim.GetCellNo());
    vector<double> x_velocity_vector(sim.GetCellNo());
    vector<double> y_velocity_vector(sim.GetCellNo());
    for(int i=0;i<sim.GetCellNo();++i){
      density_vector[i] = sim.GetCell(i).Density;
      pressure_vector[i] = sim.GetCell(i).Pressure;
      x_velocity_vector[i] = sim.GetCell(i).Velocity.x;
      y_velocity_vector[i] = sim.GetCell(i).Velocity.y;
    }
    write_std_vector_to_hdf5(file, density_vector, "density");
    write_std_vector_to_hdf5(file, pressure_vector, "pressure");
    write_std_vector_to_hdf5(file, x_velocity_vector, "x_velocity");
    write_std_vector_to_hdf5(file, y_velocity_vector, "y_velocity");
  }
}
