#include "../../misc/hdf5_utils.hpp"
#include "hdf5_diagnostics_3d.hpp"
#include <functional>

using namespace H5;
using std::function;

void write_snapshot_to_hdf5
(const HDSim3D& sim,
 const string& fname)
{
  H5File file(H5std_string(fname), H5F_ACC_TRUNC);
  Group geometry = file.createGroup("/geometry");
  Group hydrodynamic = file.createGroup("/hydrodynamic");
  
  // General
  write_std_vector_to_hdf5
    (file,
     vector<double>(1, sim.getTime()),
     "time");
  write_std_vector_to_hdf5
    (file,
     vector<double>(1, static_cast<double>(sim.getCycle())),
     "cycle");

  // Geometry
  {
    typedef function<double(const Vector3D&)> my_func;
    vector<pair<string, my_func> > fields = 
      {{"x_coordinate", [](const Vector3D& r){return r.x;}},
       {"y_coordinate", [](const Vector3D& r){return r.y;}},
       {"z_coordinate", [](const Vector3D& r){return r.z;}}};
    for_each
      (fields.begin(), 
       fields.end(),
       [&sim, &geometry]
       (const pair<string, my_func >& f)
       {
	 vector<double> buffer(sim.getTesselation().getMeshPoints().size());
	 transform
	   (sim.getTesselation().getMeshPoints().begin(),
	    sim.getTesselation().getMeshPoints().end(),
	    buffer.begin(),
	    f.second);
	 write_std_vector_to_hdf5
	   (geometry,
	    buffer,
	    f.first);
       });
  }

  // Hydrodynamics
  {
    typedef function<double(const ComputationalCell3D&)> my_func;
    vector<pair<string, my_func> >
      fields
      ({{"density", 
	 [](const ComputationalCell3D& c){return c.density;}},
	{"pressure", 
	 [](const ComputationalCell3D& c){return c.pressure;}},
	{"x_velocity", 
	 [](const ComputationalCell3D& c){return c.velocity.x;}},
	{"y_velocity", 
	 [](const ComputationalCell3D& c){return c.velocity.y;}},
	{"z_velocity",
	 [](const ComputationalCell3D& c){return c.velocity.z;}}});
    for_each
      (fields.begin(),
       fields.end(),
       [&sim, &hydrodynamic](const pair<string, my_func>& f)
       {
	 vector<double> buffer(sim.getCells().size());
	 transform
	   (sim.getCells().begin(),
	    sim.getCells().end(),
	    buffer.begin(),
	    f.second);
	 write_std_vector_to_hdf5
	   (hydrodynamic,
	    buffer,
	    f.first);
       });
  }
}
