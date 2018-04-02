#include "LagrangianExtensiveUpdater.hpp"
#include "../../misc/utils.hpp"
#include <iostream>
LagrangianExtensiveUpdater::LagrangianExtensiveUpdater(LagrangianFlux const& fc, ExtensiveUpdater const& beu)
	:fc_(fc),beu_(beu) {}

void LagrangianExtensiveUpdater::operator()
(const vector<Extensive>& fluxes,
const PhysicalGeometry& pg,
const Tessellation& tess,
const double dt,
const CacheData& cd,
const vector<ComputationalCell>& cells,
vector<Extensive>& extensives,
double time, TracerStickerNames const& ts) const
{ 
	beu_(fluxes, pg, tess, dt, cd, cells, extensives,time, ts);
	const vector<Edge>& edge_list = tess.getAllEdges();
	fc_.edge_vel_.resize(edge_list.size(),0);
	fc_.ws_.resize(edge_list.size(),0);
	size_t N = static_cast<size_t>(tess.GetPointNo());
	vector<double> dV(N, 0);
	vector<double> dv_ws(dV);
	
	vector<double> areas(edge_list.size(),0),
	  volumes(static_cast<size_t>(tess.GetPointNo()),0);
	for (size_t i = 0; i<edge_list.size(); ++i)
		areas[i]=cd.areas[i];
	for (size_t i = 0; i<static_cast<size_t>(tess.GetPointNo()); ++i)
		volumes[i]=cd.volumes[i];
	
	for (size_t i = 0; i<edge_list.size(); ++i)
	{
		const Edge& edge = edge_list[i];
		if (edge.neighbors.first < tess.GetPointNo())
		{
			dV[static_cast<size_t>(edge.neighbors.first)] += cd.areas[i] * dt*fc_.edge_vel_[i];
			dv_ws[static_cast<size_t>(edge.neighbors.first)] += cd.areas[i] * dt*fc_.ws_[i];
		}
		if (edge.neighbors.second < tess.GetPointNo())
		{
			dV[static_cast<size_t>(edge.neighbors.second)] -= cd.areas[i] * dt*fc_.edge_vel_[i];
			dv_ws[static_cast<size_t>(edge.neighbors.second)] -= cd.areas[i] * dt*fc_.ws_[i];
		}
	}
	for(size_t i=0;i<N;++i)
	{
		//if(cells[i].tracers[t_index]<2.5e8)
		{
			assert((cd.volumes[i]+dV[i]+dv_ws[i])>0);
			double a = (cd.volumes[i]+dV[i])/(cd.volumes[i]+dV[i]+dv_ws[i]);
			assert(std::isfinite(a));
			assert(a>0);
			extensives[i]*=a;
		}
	}
//	if(rank==0)
	//	std::cout<<"finished beu3"<<std::endl;
	/*
	for (size_t i = 0; i < edge_list.size(); ++i)
	{
		double ws = fc_.ws_[i];
		Edge const& edge = tess.GetEdge(static_cast<int>(i));
		if (tess.GetOriginalIndex(edge.neighbors.first) == tess.GetOriginalIndex(edge.neighbors.second))
			continue;
		double L = cd.areas[i];
		if (edge.neighbors.first < tess.GetPointNo() && cells[edge.neighbors.first].tracers[t_index]<2.75e8)
		{
			ReplaceExtensive(toadd, extensives_local[static_cast<size_t>(edge.neighbors.first)]);
			toadd *= (L*(ws)*dt / (cd.volumes[edge.neighbors.first] + dV[edge.neighbors.first] + dv_ws[edge.neighbors.first]));
			extensives[static_cast<size_t>(edge.neighbors.first)] += toadd; 
		}
		if (edge.neighbors.second < tess.GetPointNo() && cells[edge.neighbors.first].tracers[t_index]<2.75e8)
		{
			ReplaceExtensive(toadd, extensives_local[static_cast<size_t>(edge.neighbors.second)]);
			toadd *= (L*(ws)*dt / (cd.volumes[edge.neighbors.second] + dV[edge.neighbors.second] + dv_ws[edge.neighbors.second]));
			extensives[static_cast<size_t>(edge.neighbors.second)] += toadd;
		}
	}*/
	fc_.Reset();
}
