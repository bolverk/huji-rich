#include "ConstantPrimitiveEvolution.hpp"

ConstantPrimitiveEvolution::ConstantPrimitiveEvolution(bool mass_count,int n):
mass_count_(mass_count),mass_flux(0),mass_fluxt(0),N_(n)
{}

ConstantPrimitiveEvolution::~ConstantPrimitiveEvolution(void)
{}

Conserved ConstantPrimitiveEvolution::CalcFlux(Tessellation const* tessellation,
	vector<Primitive> const& cells,	double dt,
	SpatialReconstruction* interpolation,Edge const& edge,
	Vector2D const& facevelocity,RiemannSolver const& rs,int index,
	HydroBoundaryConditions const* bc,double time,vector<vector<double> > const& tracers)
{
	if(bc->IsBoundary(edge,tessellation))
		return bc->CalcFlux(tessellation,cells,facevelocity,edge,interpolation,dt,time);
	else
	{
		Vector2D normaldir = tessellation->GetMeshPoint(edge.GetNeighbor(1))-
			tessellation->GetMeshPoint(edge.GetNeighbor(0));

		Vector2D paraldir = edge.GetVertex(1) - edge.GetVertex(0);

		Primitive left = interpolation->Interpolate
			(tessellation, cells, dt, edge, 0,InBulk);
		Primitive right = interpolation->Interpolate
			(tessellation, cells, dt, edge, 1,InBulk);
		Conserved res(FluxInBulk(normaldir,paraldir,left,right,facevelocity,rs));
		int n0=edge.GetNeighbor(0);
		int n1=edge.GetNeighbor(1);
		// Do not allow outflow from this region
		/*if(((n0<N_&&res.Mass>0)||n1<N_&&res.Mass<0)&&(n0>=N_||n1>=N_))
		{
			return Conserved();
		}*/
		if(mass_count_)
		{
			if((n0>=N_&&tracers[n0][1]>1e-6)||(n1>=N_&&tracers[n1][1]>1e-6))
			{
				if(n0==index)
				{
					mass_flux+=res.Mass*edge.GetLength()*dt;
					mass_fluxt+=res.Mass*edge.GetLength()*dt*tracers[n1][1];
				}
				else
				{
					mass_flux-=res.Mass*edge.GetLength()*dt;
					mass_fluxt-=res.Mass*edge.GetLength()*dt*tracers[n0][1];
				}
			}
		}		
		return res;
	}
}

Primitive ConstantPrimitiveEvolution::UpdatePrimitive
(vector<Conserved> const& /*conservedintensive*/,
 EquationOfState const* /*eos*/,
 vector<Primitive>& cells,int index)
{
	Primitive res=cells[index];
	return res;
}

/*void ConstantPrimitiveEvolution::ReadMassFlux(string loc)const
{
	FILE* pFile;
    pFile = fopen(loc.c_str(), "rb");
	fread((void*)&mass_flux,sizeof(double),1,pFile);
	fread((void*)&mass_fluxt,sizeof(double),1,pFile);
	fclose(pFile);
*/
	/*
	fstream myFile(loc.c_str(),ios::in | ios::binary);
	myFile.read((char*)&mass_flux,sizeof(double));
	myFile.read((char*)&mass_fluxt,sizeof(double));
	cout<<"Read mass flux as "<<mass_flux<<endl;
	cout<<"Read tracer mass flux as "<<mass_fluxt<<endl;
	myFile.close();
	*/
//}

/*void ConstantPrimitiveEvolution::WriteMassFlux(string loc)const
  {*/
	/*
	fstream myFile(loc.c_str(),ios::out | ios::binary);
	myFile.write((char*)&mass_flux,sizeof(double));
	myFile.write((char*)&mass_fluxt,sizeof(double));
	myFile.close();
	*/
	  /*	FILE* pFile;
    pFile = fopen(loc.c_str(), "wb");
	fwrite((void*)&mass_flux,sizeof(double),1,pFile);
	fwrite((void*)&mass_fluxt,sizeof(double),1,pFile);
	fclose(pFile);
	}*/

double ConstantPrimitiveEvolution::GetMassFlux(void) const
{
	return mass_flux;
}
