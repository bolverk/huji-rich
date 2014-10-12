#include "mpi_macro.hpp"
#ifdef RICH_MPI

#include "marshal.hpp"

namespace
{
	void SendRecvHydro(vector<int> const& sentprocs,
		vector<vector<int> > const&
		sentcells,vector<Primitive> const& cells,
		vector<vector<double> > const& tracers,
		bool traceractive,
		EquationOfState const& eos,
		vector<vector<Primitive> > &padd,
		vector<vector<vector<double> > > &tadd)
	{
		int nlist=(int)sentprocs.size();
		const int rank = get_mpi_rank();
		const int worldsize = get_mpi_size();
		const vector<int> procorder=GetProcOrder(rank,worldsize);
		int n=worldsize-1;
		padd.resize(nlist);
		tadd.resize(nlist);
		// Send the data
		for(int i=0;i<n;++i)
		{
			int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])
				-sentprocs.begin();
			if(index<nlist)
			{
				if(rank<procorder[i])
				{
					vector<Primitive> ptemp=VectorValues(cells,sentcells[index]);
					MPI_SendVectorPrimitive(ptemp,procorder[i],0,MPI_COMM_WORLD);
					MPI_RecvVectorPrimitive(padd[index],procorder[i],0,MPI_COMM_WORLD,eos);
					if(traceractive)
					{
						vector<vector<double> > ttemp=VectorValues(tracers,sentcells[index]);
						MPI_SendVectorTracer(ttemp,procorder[i],0,MPI_COMM_WORLD);
						MPI_RecvVectorTracer(tadd[index],procorder[i],0,MPI_COMM_WORLD,
							(int)tracers[0].size());
					}
				}
				else
				{
					MPI_RecvVectorPrimitive(padd[index],procorder[i],0,MPI_COMM_WORLD,eos);
					vector<Primitive> ptemp=VectorValues(cells,sentcells[index]);
					MPI_SendVectorPrimitive(ptemp,procorder[i],0,MPI_COMM_WORLD);
					if(traceractive)
					{
						MPI_RecvVectorTracer(tadd[index],procorder[i],0,MPI_COMM_WORLD,
							(int)tracers[0].size());
						vector<vector<double> > ttemp=VectorValues(tracers,sentcells[index]);
						MPI_SendVectorTracer(ttemp,procorder[i],0,MPI_COMM_WORLD);
					}
				}
			}
		}
	}

	void GradVectorToDouble(vector<ReducedPrimitiveGradient2D> const& vec,
		vector<double> &res)
	{
		int n=(int)vec.size();
		int gradlength=8;
		if(!vec[0].tracers.empty())
			gradlength+=vec[0].tracers.size();
		if(n*gradlength!=(int)res.size())
		{
			UniversalError eo("Sizes do not match in GradVectorToDouble");
			throw eo;
		}
		for(int i=0;i<n;++i)
		{
			res[gradlength*i]=vec[i].density.x;
			res[gradlength*i+1]=vec[i].density.y;
			res[gradlength*i+2]=vec[i].pressure.x;
			res[gradlength*i+3]=vec[i].pressure.y;
			res[gradlength*i+4]=vec[i].xvelocity.x;
			res[gradlength*i+5]=vec[i].xvelocity.y;
			res[gradlength*i+6]=vec[i].yvelocity.x;
			res[gradlength*i+7]=vec[i].yvelocity.y;
			for(int j=0;j<(gradlength-8)/2;++j)
			{
				res[gradlength*i+j*2+8]=vec[i].tracers[j].x;
				res[gradlength*i+j*2+9]=vec[i].tracers[j].y;
			}
		}
	}

	void DoubleVectorToGrad(vector<ReducedPrimitiveGradient2D> &vec,
		vector<double> const& temp,int gradlength)
	{
		int n=(int)temp.size()/gradlength;
		if(n!=(int)vec.size())
		{
			UniversalError eo("Sizes do not match in DoubleVectorToGrad");
			throw eo;
		}
		for(int i=0;i<n;++i)
		{
			vec[i].density.x=temp[gradlength*i];
			vec[i].density.y=temp[gradlength*i+1];
			vec[i].pressure.x=temp[gradlength*i+2];
			vec[i].pressure.y=temp[gradlength*i+3];
			vec[i].xvelocity.x=temp[gradlength*i+4];
			vec[i].xvelocity.y=temp[gradlength*i+5];
			vec[i].yvelocity.x=temp[gradlength*i+6];
			vec[i].yvelocity.y=temp[gradlength*i+7];
			for(int j=0;j<(gradlength-8)/2;++j)
			{
				Vector2D vtemp(temp[gradlength*i+2*j+8],temp[gradlength*i+2*j+9]);
				vec[i].tracers.push_back(vtemp);
			}
		}
	}

	void DoubleVectorToTracer(vector<vector<double> > &tracer,vector<double>
		const& data,int tracerlength)
	{
		int n=(int)data.size()/tracerlength;
		if(n!=(int)tracer.size())
		{
			UniversalError eo("Sizes do not match in DoubleVectorToTracer");
			throw eo;
		}
		for(int i=0;i<n;++i)
		{
			tracer[i].resize(tracerlength);
			for(int j=0;j<tracerlength;++j)
				tracer[i][j]=data[i*tracerlength+j];
		}
	}

	void TracerVectorToDouble(vector<vector<double> > const& ttemp,vector<double>
		&dtracer)
	{
		int tracerlength=(int)ttemp[0].size();
		int n=(int)ttemp.size();
		if(n*tracerlength!=(int)dtracer.size())
		{
			UniversalError eo("Sizes do not match in TracerVectorToDouble");
			throw eo;
		}
		for(int i=0;i<n;++i)
			for(int j=0;j<tracerlength;++j)
				dtracer[i*tracerlength+j]=ttemp[i][j];
	}

	void DoubleVectorToPrimitve(vector<double> const& temp,vector<Primitive> &vec,
		EquationOfState const& eos)
	{
		int ntotal=(int) (temp.size()/4);
		if(ntotal!=(int)vec.size())
		{
			UniversalError eo("Sizes do not match in DoubleVectorToPrimitve");
			throw eo;
		}
		for(int i=0;i<ntotal;++i)
		{
			Primitive ptemp(temp[4*i],temp[4*i+1],Vector2D(temp[4*i+2],temp[4*i+3]),0,0);
			ptemp.Energy=eos.dp2e(ptemp.Density,ptemp.Pressure);
			ptemp.SoundSpeed=eos.dp2c(ptemp.Density,ptemp.Pressure);
			vec[i]=ptemp;
		}
	}

	void PrimitiveVectorToDouble(vector<Primitive> const& vec,vector<double> &res)
	{
		int n=(int)vec.size();
		if(n*4!=(int)res.size())
		{
			UniversalError eo("Sizes do not match in PrimitiveVectorToDouble");
			throw eo;
		}
		for(int i=0;i<n;++i)
		{
			res[4*i]=vec[i].Density;
			res[4*i+1]=vec[i].Pressure;
			res[4*i+2]=vec[i].Velocity.x;
			res[4*i+3]=vec[i].Velocity.y;
		}
	}

	void DoubleVectorToVector2D(vector<Vector2D> &vec,vector<double> const&
		data)
	{
		int n=(int)data.size()/2;
		if(n!=(int)vec.size())
		{
			UniversalError eo("Sizes do not match in DoubleVectorToVector2D");
			throw eo;
		}
		for(int i=0;i<n;++i)
			vec[i]=Vector2D(data[2*i],data[2*i+1]);
	}

	void Vector2DVectorToDouble(vector<double> &res,vector<Vector2D> const& data)
	{
		int n=(int)data.size();
		if(n*2!=(int)res.size())
		{
			UniversalError eo("Sizes do not match in Vector2DVectorToDouble");
			throw eo;
		}
		for(int i=0;i<n;++i)
		{
			res[2*i]=data[i].x;
			res[2*i+1]=data[i].y;
		}
	}

	int FindLoc(vector<int> const& vec,int data,int occur)
	{
		for(int i=0;i<(int)vec.size();++i)
		{
			if(vec[i]==data)
				if(occur==0)
					return i;
				else
					--occur;
		}
		UniversalError eo("Couldn't find number in vector in mpi_macro::FindLoc");
		eo.AddEntry("data to find",data);
		eo.AddEntry("Occurance number",occur);
		throw eo;
	}
}

int GetTotalPointNumber(Tessellation const& tess)
{
	int myN=tess.GetPointNo();
	int res;
	MPI_Allreduce(&myN,&res,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	return res;
}

int get_mpi_rank(void)
{
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
}

int get_mpi_size(void)
{
	int size = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	return size;
}

int MPI_Send_Vector2D(Vector2D const& vec,int dest,int tag, MPI_Comm comm)
{
	double temp[2];
	temp[0]=vec.x;
	temp[1]=vec.y;
	int err=MPI_Send(&temp[0],2,MPI_DOUBLE,dest,tag,comm);
	return err;
}

int MPI_VectorSend_Vector2D(vector<Vector2D> const& vec,int dest,int tag, MPI_Comm comm)
{
	int n=(int)vec.size();
	if(vec.empty())
	{
		double temp=0;
		int err=MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
		return err;
	}
	vector<double> temp(n*2);
	for(int i=0;i<n;++i)
	{
		temp[2*i]=vec[i].x;
		temp[2*i+1]=vec[i].y;
	}
	int err;
	err=MPI_Send(&temp[0],2*n,MPI_DOUBLE,dest,tag,comm);
	return err;
}

int MPI_Recv_Vector2D(Vector2D &vec,int source, int tag, MPI_Comm comm)
{
	double temp[2];
	int err=0;
	err=MPI_Recv(&temp[0],2,MPI_DOUBLE,source,tag,comm,MPI_STATUS_IGNORE);
	vec.x=temp[0];
	vec.y=temp[1];
	return err;
}

int MPI_VectorRecv_Vector2D(vector<Vector2D> &vec,int source, int tag, MPI_Comm comm)
{
	vec.clear();
	MPI_Status status;
	MPI_Probe(source,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	if(status.MPI_TAG==1)
	{
		double temp;
		int err=MPI_Recv(&temp,1,MPI_DOUBLE,source,1,comm,&status);
		return err;
	}
	int n=0;
	MPI_Get_count(&status,MPI_DOUBLE,&n);
	vector<double> temp(n);
	vec.resize(n/2);
	n/=2;
	int err=0;
	err=MPI_Recv(&temp[0],2*n,MPI_DOUBLE,source,tag,comm,MPI_STATUS_IGNORE);
	for(int i=0;i<n;++i)
	{
		vec[i].x=temp[2*i];
		vec[i].y=temp[2*i+1];
	}
	return err;
}

int MPI_VectorBcast_Vector2D(vector<Vector2D> &vec,int root, MPI_Comm comm,int rank)
{
	int n=(int)vec.size();
	vector<double> temp(n*2);
	if(rank==root)
	{
		for(int i=0;i<n;++i)
		{
			temp[2*i]=vec[i].x;
			temp[2*i+1]=vec[i].y;
		}
	}
	int err=MPI_Bcast(&temp[0],n*2,MPI_DOUBLE,root,comm);
	if(rank!=root)
	{
		for(int i=0;i<n;++i)
		{
			vec[i].x=temp[2*i];
			vec[i].y=temp[2*i+1];
		}
	}
	return err;
}

bool PointInsideCell(Tessellation const& tess,int cell_index,Vector2D const & point)
{
	boost::array<Vector2D,3> triangle;
	vector<int>const& cell_edges=tess.GetCellEdges(cell_index);
	triangle[0]=tess.GetEdge(cell_edges[0]).vertices.first;
	triangle[1]=tess.GetEdge(cell_edges[0]).vertices.second;
	triangle[2]=tess.GetCellCM(cell_index);
	double temp=orient2d(triangle);
	triangle[2]=point;
	return (orient2d(triangle)*temp)>0;
}

/*void BuildTree(ANNkd_tree *&tree,ANNpointArray &treePoints,
Tessellation const& tess)
{
int N=tess.GetPointNo();
treePoints=annAllocPts(N,2);
for(int i=0;i<N;++i)
{
treePoints[i][0]=tess.GetMeshPoint(i).x;
treePoints[i][1]=tess.GetMeshPoint(i).y;
}
tree=new ANNkd_tree(treePoints,N,2,1,ANN_KD_SUGGEST);
return;
}

void KillTree(ANNkd_tree *&tree,ANNpointArray &treePoints)
{
annDeallocPts(treePoints);
delete tree;
annClose();
}

int FindContainingCell(ANNkd_tree *tree,Vector2D const& point)
{
ANNpoint queryPt;
queryPt=annAllocPt(2);
queryPt[0]=point.x;
queryPt[1]=point.y;
ANNidxArray nnIdx; // near neighbor indices
ANNdistArray dists; // near neighbor distances
nnIdx = new ANNidx[1]; // allocate near neigh indices
dists = new ANNdist[1]; // allocate near neighbor dists
tree->annkSearch(queryPt,1,nnIdx,dists);
int res=nnIdx[0];
delete [] nnIdx;
delete [] dists;
annDeallocPt(queryPt);
return res;
}
*/
void ConvertVector2DToDouble(vector<Vector2D> const& vec,vector<double> &res)
{
	int n=(int) vec.size();
	res.resize(n*2);
	for(int i=0;i<n;++i)
	{
		res[2*i]=vec[i].x;
		res[2*i+1]=vec[i].y;
	}
}

vector<Vector2D> MPI_MassSendRecvVectorVector2D
	(vector<vector<Vector2D> > const& tosend,
	vector<int> const& proclist,vector<int> const& procorder)
{
	vector<Vector2D> res;
	int n=(int)procorder.size();
	int nlist=(int)proclist.size();
	const int rank = get_mpi_rank();
	for(int i=0;i<n;++i)
	{
		int index=find(proclist.begin(),proclist.end(),procorder[i])-proclist.begin();
		// Do we talk with this processor?
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				if(tosend[index].empty())
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD);
				}
				else
					MPI_VectorSend_Vector2D(tosend[index],procorder[i],0,
					MPI_COMM_WORLD);
				MPI_Status status;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				if(status.MPI_TAG==0)
				{
					vector<Vector2D> tempres;
					MPI_VectorRecv_Vector2D(tempres,procorder[i],0,MPI_COMM_WORLD);
					res.insert(res.end(),tempres.begin(),tempres.end());
				}
				else
				{
					double temp;
					MPI_Recv(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,&status);
				}
			}
			else
			{
				MPI_Status status;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				if(status.MPI_TAG==0)
				{
					vector<Vector2D> tempres;
					MPI_VectorRecv_Vector2D(tempres,procorder[i],0,MPI_COMM_WORLD);
					res.insert(res.end(),tempres.begin(),tempres.end());
				}
				else
				{
					double temp;
					MPI_Recv(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,&status);
				}
				if(tosend[index].empty())
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD);
				}
				else
					MPI_VectorSend_Vector2D(tosend[index],procorder[i],0,
					MPI_COMM_WORLD);
			}
		}
	}
	return res;
}

vector<int> GetProcOrder(int rank,int worldsize)
{
	if(worldsize==1)
		return vector<int> ();
	vector<int> procorder(worldsize-1);
	if(rank==0)
	{
		for(int i=0;i<(worldsize-1);++i)
			procorder[i]=i+1;
	}
	if(rank!=(worldsize-1))
	{
		for(int i=0;i<(worldsize-1);++i)
		{
			int temp=(i-rank+worldsize)%(worldsize-1);
			if(temp!=rank)
				procorder[i]=temp;
			else
				procorder[i]=worldsize-1;
		}
	}
	else
	{
		if(worldsize>1)
		{
			int half=(int)(worldsize/2);
			for(int i=0;i<(worldsize-1);++i)
			{
				if(i%2==0)
					procorder[i]=half+i/2;
				else
					procorder[i]=(int)(i/2)+1;
			}
			procorder[worldsize-2]=0;
		}
	}
	return procorder;
}

namespace {
	template<class T, class S> vector<T> mass_static_cast(const vector<S>& source)
	{
		vector<T> res(source.size());
		for(size_t i=0;i<source.size();++i)
			res[i] = static_cast<T>(source[i]);
		return res;
	}

	class ExtensiveCommunicator: public Communication
	{
	public:

		ExtensiveCommunicator(const vector<Conserved>& c_to_send,
			const vector<vector<double> >& t_to_send,
			const vector<size_t>& cei_to_send,int ntracer):
		c_to_send_(c_to_send), t_to_send_(t_to_send),
			cei_to_send_(cei_to_send), c_received_(),
			t_received_(), cei_received_(),ntracer_(ntracer) {}

		void sendInfo(int address)
		{
			MPI_SendVectorConserved(c_to_send_,address,0,MPI_COMM_WORLD);
			MPI_SendVectorTracer(t_to_send_,address,0,MPI_COMM_WORLD);
			if(!cei_to_send_.empty())
			{
				vector<unsigned> buf = mass_static_cast<unsigned,size_t>(cei_to_send_);
				MPI_Send(&buf[0],buf.size(),MPI_UNSIGNED,address,0,MPI_COMM_WORLD);
			}
			else
			{
				int temp = 0;
				MPI_Send(&temp,1,MPI_UNSIGNED,address,1,MPI_COMM_WORLD);
			}
		}

		void recvInfo(int address)
		{
			MPI_RecvVectorConserved(c_received_,address,0,MPI_COMM_WORLD);
			MPI_RecvVectorTracer(t_received_,address,0,MPI_COMM_WORLD,ntracer_);
			MPI_Status status;
			MPI_Probe(address,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			if(status.MPI_TAG==1)
			{
				int temp;
				MPI_Recv(&temp,1,MPI_UNSIGNED,address,1,MPI_COMM_WORLD,&status);
			}
			else
			{
				int count;
				MPI_Get_count(&status,MPI_UNSIGNED,&count);
				vector<unsigned> buf(count);
				MPI_Recv(&buf[0],count,MPI_UNSIGNED,address,0,MPI_COMM_WORLD,&status);
				cei_received_ = mass_static_cast<size_t,unsigned>(buf);
			}
		}

		const vector<Conserved>& getReceivedConserved(void) const
		{
			return c_received_;
		}

		const vector<vector<double> >& getReceivedTracers(void) const
		{
			return t_received_;
		}

		const vector<size_t>& getReceivedCustomEvolutionIndices(void) const
		{
			return cei_received_;
		}

	private:
		const vector<Conserved> c_to_send_;
		const vector<vector<double> > t_to_send_;
		const vector<size_t> cei_to_send_;
		const int ntracer_;
		vector<Conserved> c_received_;
		vector<vector<double> > t_received_;
		vector<size_t> cei_received_;
	};
}

void SendRecvExtensive(vector<Conserved> const& cons,vector<vector<double> > const&
	tracers,vector<size_t> const& customevolutions,vector<vector<int> > const& sentcells,
	vector<int> const& sentprocs,vector<Conserved> &ptoadd,vector<vector<double> >
	&ttoadd,vector<size_t> &ctoadd)
{
	bool traceractive=!tracers.empty();
	int ntracer=0;
	if(traceractive)
		traceractive=!tracers[0].empty();
	if(traceractive)
		ntracer=(int)tracers[0].size();

	// Take care of self send hydro
	int nlist=(int)sentprocs.size();
	ptoadd.clear();
	ttoadd.clear();
	ctoadd.clear();

	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);

	int n=worldsize-1;
	for(int i=0;i<n;++i)
	{
		int index=find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
		// Do we talk with this processor?
		if(index<nlist)
		{
			ExtensiveCommunicator extensive_communicator(VectorValues(cons,sentcells[index]),
				VectorValues(tracers,sentcells[index]),
				VectorValues(customevolutions,sentcells[index]),ntracer);
			marshal_communication(extensive_communicator,procorder[i],rank<procorder[i]);
			insert_all_to_back(ptoadd,extensive_communicator.getReceivedConserved());
			insert_all_to_back(ttoadd,extensive_communicator.getReceivedTracers());
			insert_all_to_back(ctoadd,extensive_communicator.getReceivedCustomEvolutionIndices());
		}
	}
}

void SendRecvShockedStatus(vector<char> const& shockedcells,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<char> &btoadd)
{
	// Take care of self send cells
	int nlist=(int)sentprocs.size();
	btoadd.clear();

	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);

	int n=worldsize-1;
	for(int i=0;i<n;++i)
	{
		int index=find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
		// Do we talk with this processor?
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<char> btemp=VectorValues(shockedcells,sentcells[index]);
				if(!btemp.empty())
					MPI_Send(&btemp[0],(int)btemp.size(),MPI_CHAR,procorder[i],0,
					MPI_COMM_WORLD);
				else
				{
					char bbtemp=0;
					MPI_Send(&bbtemp,1,MPI_CHAR,procorder[i],1,
						MPI_COMM_WORLD);
				}
				MPI_Status status;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				if(status.MPI_TAG==1)
				{
					char bbtemp;
					MPI_Recv(&bbtemp,1,MPI_CHAR,procorder[i],1,MPI_COMM_WORLD,
						&status);
				}
				else
				{
					int count;
					MPI_Get_count(&status,MPI_CHAR,&count);
					vector<char> toadd(count);
					MPI_Recv(&toadd[0],count,MPI_CHAR,procorder[i],0,
						MPI_COMM_WORLD,&status);
					btoadd.insert(btoadd.end(),toadd.begin(),toadd.end());
				}
			}
			else
			{
				MPI_Status status;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				if(status.MPI_TAG==1)
				{
					char bbtemp;
					MPI_Recv(&bbtemp,1,MPI_CHAR,procorder[i],1,MPI_COMM_WORLD,
						&status);
				}
				else
				{
					int count;
					MPI_Get_count(&status,MPI_CHAR,&count);
					vector<char> toadd(count);
					MPI_Recv(&toadd[0],count,MPI_CHAR,procorder[i],0,
						MPI_COMM_WORLD,&status);
					btoadd.insert(btoadd.end(),toadd.begin(),toadd.end());
				}
				vector<char> btemp=VectorValues(shockedcells,sentcells[index]);
				if(!btemp.empty())
					MPI_Send(&btemp[0],(int)btemp.size(),MPI_CHAR,procorder[i],0,
					MPI_COMM_WORLD);
				else
				{
					char bbtemp=0;
					MPI_Send(&bbtemp,1,MPI_CHAR,procorder[i],1,
						MPI_COMM_WORLD);
				}
			}
		}
	}
}

void SendRecvVectorDouble(vector<double> const& vec,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<double> &toadd)
{
	// Take care of self send cells
	int nlist=(int)sentprocs.size();
	toadd.clear();

	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);

	int n=worldsize-1;
	for(int i=0;i<n;++i)
	{
		int index=find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
		// Do we talk with this processor?
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<double> temp=VectorValues(vec,sentcells[index]);
				if(!temp.empty())
					MPI_Send(&temp[0],(int)temp.size(),MPI_DOUBLE,procorder[i],0,
					MPI_COMM_WORLD);
				else
				{
					double bbtemp=0;
					MPI_Send(&bbtemp,1,MPI_DOUBLE,procorder[i],1,
						MPI_COMM_WORLD);
				}
				MPI_Status status;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				if(status.MPI_TAG==1)
				{
					double bbtemp;
					MPI_Recv(&bbtemp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,
						&status);
				}
				else
				{
					int count;
					MPI_Get_count(&status,MPI_DOUBLE,&count);
					vector<double> ttoadd(count);
					MPI_Recv(&ttoadd[0],count,MPI_DOUBLE,procorder[i],0,
						MPI_COMM_WORLD,&status);
					toadd.insert(toadd.end(),ttoadd.begin(),ttoadd.end());
				}
			}
			else
			{
				MPI_Status status;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				if(status.MPI_TAG==1)
				{
					double bbtemp;
					MPI_Recv(&bbtemp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD,
						&status);
				}
				else
				{
					int count;
					MPI_Get_count(&status,MPI_DOUBLE,&count);
					vector<double> ttoadd(count);
					MPI_Recv(&ttoadd[0],count,MPI_DOUBLE,procorder[i],0,
						MPI_COMM_WORLD,&status);
					toadd.insert(toadd.end(),ttoadd.begin(),ttoadd.end());
				}
				vector<double> temp=VectorValues(vec,sentcells[index]);
				if(!temp.empty())
					MPI_Send(&temp[0],(int)temp.size(),MPI_DOUBLE,procorder[i],0,
					MPI_COMM_WORLD);
				else
				{
					double bbtemp=0;
					MPI_Send(&bbtemp,1,MPI_DOUBLE,procorder[i],1,
						MPI_COMM_WORLD);
				}
			}
		}
	}
}

void SendRecvOldVector2D(vector<Vector2D> const& points,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<Vector2D> &toadd)
{
	int nlist=(int)sentprocs.size();
	toadd.clear();

	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);

	int n=worldsize-1;
	for(int i=0;i<n;++i)
	{
		int index=find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin();
		// Do we talk with this processor?
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<Vector2D> vtemp=VectorValues(points,sentcells[index]);
				MPI_VectorSend_Vector2D(vtemp,procorder[i],0,MPI_COMM_WORLD);
				MPI_VectorRecv_Vector2D(vtemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD);
				if(!vtemp.empty())
				{
					toadd.insert(toadd.end(),vtemp.begin(),vtemp.end());
				}
			}
			else
			{
				vector<Vector2D> vtemp;
				MPI_VectorRecv_Vector2D(vtemp,procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD);
				if(!vtemp.empty())
				{
					toadd.insert(toadd.end(),vtemp.begin(),vtemp.end());
				}
				vtemp=VectorValues(points,sentcells[index]);
				MPI_VectorSend_Vector2D(vtemp,procorder[i],0,MPI_COMM_WORLD);
			}
		}
	}
}

void KeepLocalPoints(vector<Conserved> &cons,vector<vector<double> > &tracers,
	vector<size_t> &customevolutions,vector<int> const& localpoints)
{
	bool traceractive=!tracers.empty();
	if(traceractive)
		traceractive=!tracers[0].empty();
	cons=VectorValues(cons,localpoints);
	if(traceractive)
		tracers=VectorValues(tracers,localpoints);
	customevolutions=VectorValues(customevolutions,localpoints);
}

namespace {
	class HydroCommunicator: public Communication
	{
	public:

		HydroCommunicator(const EquationOfState& eos,
			const vector<Primitive>& p_to_send,
			const vector<size_t>& cei_to_send):
		eos_(eos), p_to_send_(p_to_send),
			cei_to_send_(cei_to_send),
			p_received_(), cei_received_() {}

		void sendInfo(int address)
		{
			MPI_SendVectorPrimitive(p_to_send_,address,0,MPI_COMM_WORLD);
			if(!cei_to_send_.empty()){
				vector<unsigned> buf = mass_static_cast<unsigned,size_t>(cei_to_send_);
				MPI_Send(&buf[0],static_cast<int>(buf.size()),MPI_UNSIGNED,address,0,MPI_COMM_WORLD);
			}
			else{
				unsigned temp = 0;
				MPI_Send(&temp,1,MPI_UNSIGNED,address,1,MPI_COMM_WORLD);
			}
		}

		void recvInfo(int address)
		{
			try
			{
				MPI_RecvVectorPrimitive(p_received_,address,0,MPI_COMM_WORLD,eos_);
			}
			catch(UniversalError &eo)
			{
				eo.AddEntry("Error in HydroCommunicator while recv from cpu",address);
				throw;
			}
			int count;
			MPI_Status stat;
			MPI_Probe(address,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
			MPI_Get_count(&stat,MPI_UNSIGNED,&count);
			if(stat.MPI_TAG==0){
				vector<unsigned> buf(count);
				MPI_Recv(&buf[0],count,MPI_UNSIGNED,address,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				cei_received_ = mass_static_cast<size_t,unsigned>(buf);
			}
			else{
				unsigned temp;
				MPI_Recv(&temp,1,MPI_UNSIGNED,address,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}

		const vector<Primitive>& getReceivedPrimitives(void) const
		{
			return p_received_;
		}

		const vector<size_t> getReceivedCustomEvolutionIndices(void) const
		{
			return cei_received_;
		}

	private:
		const EquationOfState& eos_;
		const vector<Primitive> p_to_send_;
		const vector<size_t> cei_to_send_;
		vector<Primitive> p_received_;
		vector<size_t> cei_received_;
	};
}

void SendRecvHydro(vector<Primitive> &cells,
	vector<size_t> &customevolutions,vector<vector<int> > const& sentcells,
	vector<int> const& sentprocs,EquationOfState const& eos,vector<vector<int> > const& Nghost,
	int totalpoints)
{
	int nlist=(int)sentprocs.size();
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int n=worldsize-1;
	vector<vector<Primitive> > padd(nlist);
	vector<vector<size_t> > cadd(nlist);

	// Send the data
	for(int i=0;i<n;++i)
	{
		int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])
			-sentprocs.begin();
		if(index==nlist)
			continue;
		HydroCommunicator hydro_communicator(eos,
			VectorValues(cells,sentcells[index]),
			VectorValues(customevolutions,sentcells[index]));
		marshal_communication(hydro_communicator,procorder[i],rank<procorder[i]);
		padd[index] = hydro_communicator.getReceivedPrimitives();
		cadd[index] = hydro_communicator.getReceivedCustomEvolutionIndices();
	}
	// ReArrange the data
	cells.resize(totalpoints);
	customevolutions.resize(totalpoints);
	for(int i=0;i<nlist;++i)
	{
		ListExchange(cells,Nghost[i],padd[i]);
		ListExchange(customevolutions,Nghost[i],cadd[i]);
	}
}

namespace {
	class TracerCommunicator: public Communication
	{
	public:

		TracerCommunicator(const vector<vector<double> >& to_send,int dim):
		  to_send_(to_send),dim_(dim),reply_() {}

		  void sendInfo(int address)
		  {
			  MPI_SendVectorTracer(to_send_,address,0,MPI_COMM_WORLD);
		  }

		  void recvInfo(int address)
		  {
			  MPI_RecvVectorTracer(reply_,address,0,MPI_COMM_WORLD,dim_);
		  }

		  const vector<vector<double> >& getReply(void) const
		  {
			  return reply_;
		  }

	private:
		const vector<vector<double> > to_send_;
		const int dim_;
		vector<vector<double> > reply_;
	};
}

void SendRecvTracers(vector<vector<double> > &tracers,
	vector<vector<int> > const& sentcells,vector<int> const& sentprocs,
	vector<vector<int> > const& Nghost,int totalpoints)
{
	int nlist=(int)sentprocs.size();
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	vector<vector<vector<double> > > tadd(nlist);
	// Send the data
	for(int i=0;i<worldsize-1;++i)
	{
		int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])
			-sentprocs.begin();
		if(index<nlist)
		{
			TracerCommunicator tc(VectorValues(tracers,sentcells[index]),
				(int)tracers[0].size());
			marshal_communication(tc,procorder[i],
				rank<procorder[i]);
			tadd[index] = tc.getReply();
		}
	}
	// ReArrange the data
	tracers.resize(totalpoints);
	for(int i=0;i<nlist;++i)
		ListExchange(tracers,Nghost[i],tadd[i]);
}

int MPI_SendVectorPrimitive(vector<Primitive> const& vec,int dest,int tag,
	MPI_Comm comm)
{
	if(vec.empty())
	{
		double temp=0;
		return MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
	}
	int n=(int)vec.size();
	vector<double> tosend(n*4);
	for(int i=0;i<n;++i)
	{
		tosend[4*i]=vec[i].Density;
		tosend[4*i+1]=vec[i].Pressure;
		tosend[4*i+2]=vec[i].Velocity.x;
		tosend[4*i+3]=vec[i].Velocity.y;
	}
	return MPI_Send(&tosend[0],4*n,MPI_DOUBLE,dest,tag,comm);
}

void MPI_RecvVectorPrimitive(vector<Primitive> &vec,int dest,int tag,
	MPI_Comm comm,EquationOfState const& eos)
{
	vec.clear();
	MPI_Status status;
	MPI_Probe(dest,MPI_ANY_TAG,comm,&status);
	if(status.MPI_TAG==1)
	{
		double temp;
		MPI_Recv(&temp,1,MPI_DOUBLE,dest,1,comm,&status);
		return;
	}
	if(status.MPI_TAG!=tag)
	{
		UniversalError eo("recveived wrong tag in MPI_RecvVectorPrimitive");
		throw eo;
	}
	int nrecv;
	MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
	vector<double> temp(nrecv);
	MPI_Recv(&temp[0],nrecv,MPI_DOUBLE,dest,tag,comm,&status);
	int ntotal=nrecv/4;
	vec.reserve(ntotal);
	for(int i=0;i<ntotal;++i)
	{
		Primitive ptemp(temp[4*i],temp[4*i+1],Vector2D(temp[4*i+2],temp[4*i+3]),0,0);
		ptemp.Energy=eos.dp2e(ptemp.Density,ptemp.Pressure);
		ptemp.SoundSpeed=eos.dp2c(ptemp.Density,ptemp.Pressure);
		vec.push_back(ptemp);
	}
}

int MPI_SendVectorTracer(vector<vector<double> > const& vec,int dest,int tag,
	MPI_Comm comm)
{
	if(vec.empty())
	{
		double temp=0;
		return MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
	}
	int n=(int)vec.size();
	int ntracer=(int)vec[0].size();
	vector<double> tosend(n*ntracer);
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<ntracer;++j)
			tosend[ntracer*i+j]=vec[i][j];
	}
	return MPI_Send(&tosend[0],ntracer*n,MPI_DOUBLE,dest,tag,comm);
}

void MPI_RecvVectorTracer(vector<vector<double> > &vec,int dest,int tag,
	MPI_Comm comm,int ntracer)
{
	vec.clear();
	int nrecv;
	MPI_Status status;
	MPI_Probe(dest,MPI_ANY_TAG,comm,&status);
	if(status.MPI_TAG==1)
	{
		double temp;
		MPI_Recv(&temp,1,MPI_DOUBLE,dest,1,comm,&status);
		return;
	}
	if(status.MPI_TAG!=tag)
	{
		UniversalError eo("recveived wrong tag in MPI_RecvVectorTracer");
		throw eo;
	}
	MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
	vector<double> temp(nrecv);
	MPI_Recv(&temp[0],nrecv,MPI_DOUBLE,dest,MPI_ANY_TAG,comm,&status);
	int length=nrecv/ntracer;
	vec.resize(length);
	for(int i=0;i<length;++i)
	{
		vec[i].resize(ntracer);
		for(int j=0;j<ntracer;++j)
			vec[i][j]=temp[ntracer*i+j];
	}
}

int MPI_SendVectorGrad(vector<ReducedPrimitiveGradient2D> const&vec,int dest,int
	tag,MPI_Comm comm)
{
	if(vec.empty())
	{
		double temp=0;
		return MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
	}
	int n=(int)vec.size();
	int gradlength=2*(int)vec[0].tracers.size()+8;
	vector<double> tosend(n*gradlength);
	for(int i=0;i<n;++i)
	{
		tosend[gradlength*i]=vec[i].density.x;
		tosend[gradlength*i+1]=vec[i].density.y;
		tosend[gradlength*i+2]=vec[i].pressure.x;
		tosend[gradlength*i+3]=vec[i].pressure.y;
		tosend[gradlength*i+4]=vec[i].xvelocity.x;
		tosend[gradlength*i+5]=vec[i].xvelocity.y;
		tosend[gradlength*i+6]=vec[i].yvelocity.x;
		tosend[gradlength*i+7]=vec[i].yvelocity.y;
		for(int j=0;j<(gradlength-8)/2;++j)
		{
			tosend[gradlength*i+j*2+8]=vec[i].tracers[j].x;
			tosend[gradlength*i+j*2+9]=vec[i].tracers[j].y;
		}
	}
	return MPI_Send(&tosend[0],gradlength*n,MPI_DOUBLE,dest,tag,comm);
}

void MPI_RecvVectorGrad(vector<ReducedPrimitiveGradient2D> &vec,int dest,int
	tag,MPI_Comm comm,int gradlength)
{
	MPI_Status status;
	MPI_Probe(dest,MPI_ANY_TAG,comm,&status);
	if(status.MPI_TAG==1)
	{
		double temp;
		MPI_Recv(&temp,1,MPI_DOUBLE,dest,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		return;
	}
	if(status.MPI_TAG!=tag)
	{
		UniversalError eo("recveived wrong tag in MPI_RecvVectorGrad");
		throw eo;
	}
	int ntotal;
	MPI_Get_count(&status,MPI_DOUBLE,&ntotal);
	vector<double> temp(ntotal);
	MPI_Recv(&temp[0],ntotal,MPI_DOUBLE,dest,tag,comm,&status);
	int n=ntotal/gradlength;
	vec.resize(n);
	for(int i=0;i<n;++i)
	{
		vec[i].density.x=temp[gradlength*i];
		vec[i].density.y=temp[gradlength*i+1];
		vec[i].pressure.x=temp[gradlength*i+2];
		vec[i].pressure.y=temp[gradlength*i+3];
		vec[i].xvelocity.x=temp[gradlength*i+4];
		vec[i].xvelocity.y=temp[gradlength*i+5];
		vec[i].yvelocity.x=temp[gradlength*i+6];
		vec[i].yvelocity.y=temp[gradlength*i+7];
		for(int j=0;j<(gradlength-8)/2;++j)
		{
			Vector2D vtemp(temp[gradlength*i+2*j+8],temp[gradlength*i+2*j+9]);
			vec[i].tracers.push_back(vtemp);
		}
	}
}

void SendRecvVelocity(vector<Vector2D> &vel,vector<vector<int> >const& sentcells,
	vector<int> sentprocs,vector<vector<int> > const& Nghost,int totalpoints)
{
	int nlist=(int)sentprocs.size();
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int n=worldsize-1;
	vector<vector<Vector2D> > tadd(nlist);
	// Send the data
	for(int i=0;i<n;++i)
	{
		int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])
			-sentprocs.begin();
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<Vector2D> ttemp=VectorValues(vel,sentcells[index]);
				MPI_VectorSend_Vector2D(ttemp,procorder[i],0,MPI_COMM_WORLD);
				MPI_VectorRecv_Vector2D(tadd[index],procorder[i],0,MPI_COMM_WORLD);
			}
			else
			{
				MPI_VectorRecv_Vector2D(tadd[index],procorder[i],0,MPI_COMM_WORLD);
				vector<Vector2D> ttemp=VectorValues(vel,sentcells[index]);
				MPI_VectorSend_Vector2D(ttemp,procorder[i],0,MPI_COMM_WORLD);
			}
		}
	}
	// ReArrange the data
	vel.resize(totalpoints);
	for(int i=0;i<nlist;++i)
		ListExchange(vel,Nghost[i],tadd[i]);
}

void SendRecvGrad(vector<ReducedPrimitiveGradient2D> &grads,
	vector<vector<int> >const& sentcells,vector<int> sentprocs,
	vector<vector<int> > const& Nghost,int totalpoints)
{
	if(grads.empty())
		return;
	int nlist=(int)sentprocs.size();
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int n=worldsize-1;
	vector<vector<ReducedPrimitiveGradient2D> > tadd(nlist);
	int gradlength=8+2*(int)grads[0].tracers.size();
	// Send the data
	for(int i=0;i<n;++i)
	{
		int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])
			-sentprocs.begin();
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<ReducedPrimitiveGradient2D> ttemp=VectorValues(grads,sentcells[index]);
				MPI_SendVectorGrad(ttemp,procorder[i],0,MPI_COMM_WORLD);
				MPI_RecvVectorGrad(tadd[index],procorder[i],0,MPI_COMM_WORLD,gradlength);
			}
			else
			{
				vector<ReducedPrimitiveGradient2D> ttemp=VectorValues(grads,sentcells[index]);
				MPI_RecvVectorGrad(tadd[index],procorder[i],0,MPI_COMM_WORLD,gradlength);
				MPI_SendVectorGrad(ttemp,procorder[i],0,MPI_COMM_WORLD);
			}
		}
	}
	// ReArrange the data
	grads.resize(totalpoints);
	for(int i=0;i<nlist;++i)
		ListExchange(grads,Nghost[i],tadd[i]);
}

int MPI_SendVectorConserved(vector<Conserved> const& vec,int dest,int tag,
	MPI_Comm comm)
{
	if(vec.empty())
	{
		double temp=0;
		return MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
	}
	int n=(int)vec.size();
	vector<double> tosend(n*4);
	for(int i=0;i<n;++i)
	{
		tosend[4*i]=vec[i].Mass;
		tosend[4*i+1]=vec[i].Energy;
		tosend[4*i+2]=vec[i].Momentum.x;
		tosend[4*i+3]=vec[i].Momentum.y;
	}
	return MPI_Send(&tosend[0],4*n,MPI_DOUBLE,dest,tag,comm);
}

void MPI_RecvVectorConserved(vector<Conserved> &vec,int dest,int tag,
	MPI_Comm comm)
{
	vec.clear();
	MPI_Status status;
	MPI_Probe(dest,MPI_ANY_TAG,comm,&status);
	if(status.MPI_TAG==1)
	{
		double temp;
		MPI_Recv(&temp,1,MPI_DOUBLE,dest,1,comm,&status);
		return;
	}
	if(status.MPI_TAG!=tag)
	{
		UniversalError eo("recveived wrong tag in MPI_RecvVectorConserved");
		throw eo;
	}

	int nrecv;
	MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
	vector<double> temp(nrecv);
	MPI_Recv(&temp[0],nrecv,MPI_DOUBLE,dest,tag,comm,&status);
	int ntotal=nrecv/4;
	vec.reserve(ntotal);
	for(int i=0;i<ntotal;++i)
	{
		Conserved ctemp(temp[4*i],Vector2D(temp[4*i+2],temp[4*i+3]),temp[4*i+1]);
		vec.push_back(ctemp);
	}
}

void SendRecvGhostIndeces(vector<vector<int> > &GhostIndeces,vector<int>
	const& BoundaryPoints,vector<vector<int> > const& SentPoints,vector<int> const&
	SentProcs)
{
	int nprocs=(int)SentProcs.size();
	int nbound=(int)BoundaryPoints.size();
	/*
	const int rank = get_mpi_rank();
	const int ws = get_mpi_size();
	*/
	vector<MPI_Status> status(nprocs);
	vector<MPI_Request> req(nprocs);

	vector<vector<int> > tosend(nprocs),torecv(nprocs);
	vector<int> sentme,flags(nprocs,0);
	int itemp=-1;

	for(int i=0;i<nprocs;++i)
	{
		// Send the data
		if(!SentPoints[i].empty())
		{
			// Find the relevant points
			vector<int> temp(SentPoints[i]);
			sort(temp.begin(),temp.end());
			for(int j=0;j<nbound;++j)
			{
				const int index=lower_bound(temp.begin(),temp.end(),BoundaryPoints[j])
					-temp.begin();
				if(index<(int)temp.size())
					tosend[i].push_back(index);
			}
			MPI_Isend(&tosend[i][0],(int)tosend.size(),MPI_INT,SentProcs[i],0,
				MPI_COMM_WORLD,&req[i]);
		}
		else
		{
			MPI_Isend(&itemp,1,MPI_INT,SentProcs[i],1,
				MPI_COMM_WORLD,&req[i]);
		}
	}
	//Recv data
	int counter=0;
	while(counter<nprocs)
	{
		for(int i=0;i<nprocs;++i)
		{
			if(!flags[i])
			{
				MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flags[i],&status[i]);
				if(flags[i])
				{
					++counter;
					if(status[i].MPI_TAG==1)
					{
						int itemp2;
						MPI_Irecv(&itemp2,1,MPI_INT,status[i].MPI_SOURCE,1,
							MPI_COMM_WORLD,&req[i]);
					}
					else
					{
						int count;
						MPI_Get_count(&status[i],MPI_INT,&count);
						torecv[i].resize(count);
						MPI_Irecv(&torecv[i][0],count,MPI_INT,status[i].MPI_SOURCE,0,
							MPI_COMM_WORLD,&req[i]);
						sentme.push_back(status[i].MPI_SOURCE);
					}
				}
			}
		}
	}
	MPI_Waitall(nprocs,&req[0],&status[0]);
	// ReArrange the data
	vector<int> occur(SentProcs);
	sort(occur.begin(),occur.end());
	occur=unique(occur);
	GhostIndeces.resize(nprocs);
	for(int i=0;i<nprocs;++i)
	{
		int loc=lower_bound(occur.begin(),occur.end(),sentme[i])-occur.begin();
		int index=FindLoc(SentProcs,sentme[i],occur[loc]);
		++occur[loc];
		GhostIndeces[index]=torecv[i];
	}
}

namespace
{
	// Nghost is the NghostIndex
	// ghost is the index in cor of the ghost removed points
	// Sentindex is the index of the proc we are dealing with
	// returns all of the local neighboring points of the removed ghost point
	vector<vector<int> > FindLocalNeighbors(vector<int> const& Nghost,vector<int>  &ghost,
		int SentIndex,Tessellation const& tess)
	{
		vector<vector<int> > const& duplicated=tess.GetDuplicatedPoints();
		vector<int> const& Sent=duplicated[SentIndex];
		vector<int> index;
		sort_index(ghost,index);
		sort(ghost.begin(),ghost.end());
		int nreal=(int)ghost.size();
		vector<vector<int> > res(ghost.size());
		for(int i=0;i<(int)Sent.size();++i)
		{ //do i need to check from other procs??????
			vector<int> neigh=tess.GetNeighbors(Sent[i]);
			for(int j=0;j<(int)neigh.size();++j)
			{
				if(binary_search(ghost.begin(),ghost.end(),neigh[j]))
				{
					int index2=lower_bound(ghost.begin(),ghost.end(),neigh[j])-
						ghost.begin();
					res[index[index2]].push_back(Sent[i]);
				}
			}
		}
		for(int i=0;i<nreal;++i)
		{
			if(res[i].empty())
			{
				UniversalError eo("empty vector in FindLocalNeighbors");
				throw eo;
			}
			sort(res[i].begin(),res[i].end());
			res[i]=unique(res[i]);
		}
		return res;
	}

	// converts indeces to real point location
	// BoundaryRemove is the list per proc what points are removed given by their indeces in the Nghost vector
	// BoundaryNeigh, for each point in boundary remove, what are the indeces in Nghost of the local neighbors
	void GetRealNeighbors(Tessellation const& tess,vector<vector<int> > &BoundaryRemove,
		vector<vector<vector<int> > > &BoundaryNeigh,vector<vector<int> > &localneigh,
		vector<vector<int> > &ghostneigh)
	{
		int nprocs=(int)BoundaryRemove.size();
		vector<vector<int> > const& Nghost=tess.GetGhostIndeces();
		localneigh.clear();
		ghostneigh.clear();
		for(int i=0;i<nprocs;++i)
		{
			if(BoundaryRemove.empty())
				continue;
			vector<int> ghosts(BoundaryRemove[i].size());
			for(int j=0;j<(int)BoundaryRemove[i].size();++j)
			{
				vector<int> itemp;
				for(vector<int>::iterator it=BoundaryNeigh[i][j].begin();it!=
					BoundaryNeigh[i][j].end();++it)
					itemp.push_back(Nghost[i][*it]);
				sort(itemp.begin(),itemp.end());
				vector<int> rtemp(1,Nghost[i][BoundaryRemove[i][j]]);
				RemoveList(itemp,rtemp);
				itemp.insert(itemp.begin(),rtemp[0]);
				ghostneigh.push_back(itemp);
				ghosts[j]=Nghost[i][BoundaryRemove[i][j]];
			}
			vector<vector<int> > temp=FindLocalNeighbors(Nghost[i],ghosts,i,tess);
			for(int j=0;j<(int)BoundaryRemove[i].size();++j)
				localneigh.push_back(temp[j]);
		}
	}
}

void GetAMRExtensive(vector<Primitive> &rescells,
	vector<vector<double> > &restracer,
	vector<Primitive> const& cells,
	vector<vector<double> >
	const& tracers,
	bool traceractive,
	vector<vector<int> > &ToSend,
	vector<int> const& proclist,
	EquationOfState const& eos,
	vector<vector<int> > const& DuplicatedPoints,
	vector<vector<int> > const& Nghost,
	vector<int> const& ToRemove)
{
	// ToSend is the index in the Nghost
	// Send to other procs the list of points
	int nlist=(int)proclist.size();
	vector<vector<int> > recv(nlist);
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int n=worldsize-1;
	int temp;
	// Send the data
	for(int i=0;i<n;++i)
	{
		int index=Find(proclist.begin(),proclist.end(),procorder[i])
			-proclist.begin();
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				if(ToSend[index].empty())
					MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
					MPI_Send(&ToSend[index][0],(int)ToSend[index].size(),MPI_INT,
					procorder[i],0,MPI_COMM_WORLD);
				MPI_Status stat;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					recv[index].resize(count);
					MPI_Recv(&recv[index][0],count,MPI_INT,procorder[i],0,
						MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
			}
			else
			{
				MPI_Status stat;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					recv[index].resize(count);
					MPI_Recv(&recv[index][0],count,MPI_INT,procorder[i],0,
						MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if(ToSend[index].empty())
					MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
					MPI_Send(&ToSend[index][0],(int)ToSend[index].size(),MPI_INT,
					procorder[i],0,MPI_COMM_WORLD);
			}
		}
	}
	vector<vector<int> > ToSend2(nlist);
	for(int i=0;i<nlist;++i)
		for(int j=0;j<(int)recv[i].size();++j)
			ToSend2[i].push_back(DuplicatedPoints[i][recv[i][j]]);
	// recv is the index in the duplicated that was recieved
	vector<vector<Primitive> > padd;
	vector<vector<vector<double> > > tadd;
	SendRecvHydro(proclist,ToSend2,cells,tracers,traceractive,eos,padd,tadd);
	// rearrange the data
	int nadd=ElementNumber(padd);
	rescells.reserve(nadd);
	if(traceractive)
		restracer.reserve(nadd);
	// convert ToSend to their nghost value and sort them
	for(int i=0;i<nlist;++i)
	{
		if(!ToSend[i].empty())
		{
			for(int j=0;j<(int)ToSend[i].size();++j)
				ToSend[i][j]=Nghost[i][ToSend[i][j]];
			vector<int> indeces;
			sort_index(ToSend[i],indeces);
			sort(ToSend[i].begin(),ToSend[i].end());
			ReArrangeVector(padd[i],indeces);
			if(traceractive)
				ReArrangeVector(tadd[i],indeces);
		}
	}
	// rearrange the data
	for(int k=0;k<(int)ToRemove.size();++k)
	{
		for(int i=0;i<nlist;++i)
		{
			if(binary_search(ToSend[i].begin(),ToSend[i].end(),ToRemove[k]))
			{
				int index2=lower_bound(ToSend[i].begin(),ToSend[i].end(),
					ToRemove[k])-ToSend[i].begin();
				rescells.push_back(padd[i][index2]);
				if(traceractive)
					restracer.push_back(tadd[i][index2]);
				break;
			}
		}
		if((int)rescells.size()!=(k+1))
		{
			UniversalError eo("couldn't find point to remove int GetAMRExtensive");
			throw eo;
		}
	}
}

// BoundaryRemove is the list per proc what points are removed given by their indeces in the sent vector
// BoundaryNeigh, for each point in boundary remove, what are the indeces in sentvector of the local neighbors
void SendRecvBoundaryRemove(vector<vector<int> > &BoundaryRemove,
	vector<vector<vector<int> > > &BoundaryNeigh,Tessellation const& tess,
	vector<vector<int> > &localNeighbors,vector<vector<int> > &ghostneigh)
{
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	vector<int> sentprocs=tess.GetDuplicatedProcs();
	int nlist=(int)sentprocs.size();
	int n=worldsize-1;
	vector<int> proclist=tess.GetDuplicatedProcs();
	vector<vector<int> > RecvPoints(proclist.size());
	vector<vector<vector<int> > > RecvNeigh(proclist.size());
	int trash,trashr;
	MPI_Status stat;
	for(int i=0;i<n;++i)
	{
		//	cout<<"rank "<<rank<<" i "<<i<<" talk with "<<procorder[i]<<endl;
		int index=Find(sentprocs.begin(),sentprocs.end(),procorder[i])
			-sentprocs.begin();
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				if(BoundaryRemove[index].empty())
					MPI_Send(&trash,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
					MPI_Send(&BoundaryRemove[index][0],(int)BoundaryRemove[index].size(),
					MPI_INT,procorder[i],0,MPI_COMM_WORLD);
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&trashr,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					RecvPoints[index].resize(count);
					MPI_Recv(&RecvPoints[index][0],count,MPI_INT,procorder[i],0,
						MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				// create the lengths of neighbors
				vector<int> senddata,lengthr,recvdata;
				if(BoundaryRemove[index].empty())
					MPI_Send(&trash,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
				{
					vector<int> lengths(BoundaryRemove[index].size());
					for(int j=0;j<(int)BoundaryNeigh[index].size();++j)
					{
						lengths[j]=BoundaryNeigh[index][j].size();
						for(int k=0;k<lengths[j];++k)
							senddata.push_back(BoundaryNeigh[index][j][k]);
					}
					MPI_Send(&lengths[0],(int)lengths.size(),MPI_INT,procorder[i],0,
						MPI_COMM_WORLD);
				}
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&trashr,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					lengthr.resize(count);
					MPI_Recv(&lengthr[0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if(!BoundaryRemove[index].empty())
					MPI_Send(&senddata[0],(int)senddata.size(),MPI_INT,procorder[i],0,MPI_COMM_WORLD);
				if(!lengthr.empty())
				{
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					recvdata.resize(count);
					MPI_Recv(&recvdata[0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					//	RecvNeigh[index].resize(lengthr.size());
				}
				int loc=0;
				// Reorganize the data
				for(int j=0;j<(int)lengthr.size();++j)
				{
					vector<int> itemp(recvdata.begin()+loc,recvdata.begin()+loc+
						lengthr[j]);
					loc+=lengthr[j];
					RecvNeigh[index].push_back(itemp);
				}
			}
			else
			{
				vector<int> senddata,lengthr,recvdata;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&trashr,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					RecvPoints[index].resize(count);
					MPI_Recv(&RecvPoints[index][0],count,MPI_INT,procorder[i],0,
						MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if(BoundaryRemove[index].empty())
					MPI_Send(&trash,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
					MPI_Send(&BoundaryRemove[index][0],(int)BoundaryRemove[index].size(),
					MPI_INT,procorder[i],0,MPI_COMM_WORLD);
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&trashr,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					lengthr.resize(count);
					MPI_Recv(&lengthr[0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if(BoundaryRemove[index].empty())
					MPI_Send(&trash,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
				{
					vector<int> lengths(BoundaryRemove[index].size());
					for(int j=0;j<(int)BoundaryNeigh[index].size();++j)
					{
						lengths[j]=BoundaryNeigh[index][j].size();
						for(int k=0;k<lengths[j];++k)
							senddata.push_back(BoundaryNeigh[index][j][k]);
					}
					MPI_Send(&lengths[0],(int)lengths.size(),MPI_INT,procorder[i],0,
						MPI_COMM_WORLD);
				}
				if(!lengthr.empty())
				{
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					recvdata.resize(count);
					MPI_Recv(&recvdata[0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					//RecvNeigh[index].resize(lengthr.size());
				}
				if(!BoundaryRemove[index].empty())
					MPI_Send(&senddata[0],(int)senddata.size(),MPI_INT,procorder[i],0,MPI_COMM_WORLD);
				int loc=0;
				// Reorganize the data
				for(int j=0;j<(int)lengthr.size();++j)
				{
					vector<int> itemp(recvdata.begin()+loc,recvdata.begin()+loc+
						lengthr[j]);
					loc+=lengthr[j];
					RecvNeigh[index].push_back(itemp);
				}
			}
		}
	}
	GetRealNeighbors(tess,RecvPoints,RecvNeigh,localNeighbors,ghostneigh);
}

vector<int> RemoveMPINeighbors(vector<int> const& toremove,vector<double> const& merit,
	Tessellation const& tess)
{
	// remove is sorted
	// Find boundary cells
	int nremove=(int)toremove.size();
	vector<int> proclist=tess.GetDuplicatedProcs();
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int nlist=(int)proclist.size();
	vector<vector<int> > bremove(nlist);
	vector<vector<double> > bmerit(nlist);
	vector<vector<int> > sentpoints=tess.GetDuplicatedPoints();
	vector<vector<int> > sort_indeces(nlist);
	for(int i=0;i<nlist;++i)
		if(!sentpoints[i].empty())
		{
			sort_index(sentpoints[i],sort_indeces[i]);
			sort(sentpoints[i].begin(),sentpoints[i].end());
		}
		for(int i=0;i<nremove;++i)
		{
			for(int j=0;j<nlist;++j)
			{
				if(!binary_search(sentpoints[j].begin(),sentpoints[j].end(),
					toremove[i]))
					continue;
				int index=lower_bound(sentpoints[j].begin(),sentpoints[j].end(),
					toremove[i])-sentpoints[j].begin();
				if(index<(int)sentpoints[j].size())
				{
					bremove[j].push_back(sort_indeces[j][index]);
					bmerit[j].push_back(merit[i]);
				}
			}
		}

		MPI_Status status;
		// Send/Recv the data
		vector<vector<int> > recvindex(nlist); // the index in the Nghost vector
		vector<vector<double> > recvmerit(nlist);
		int temp;
		for(int i=0;i<(int)procorder.size();++i)
		{
			int index=Find(proclist.begin(),proclist.end(),procorder[i])
				-proclist.begin();
			if(index<nlist)
			{
				if(rank<procorder[i])
				{
					if(bremove[index].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
						MPI_Send(&bremove[index][0],(int)bremove[index].size(),
						MPI_INT,procorder[i],0,MPI_COMM_WORLD);
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					if(status.MPI_TAG==1)
						MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					else
					{
						int count;
						MPI_Get_count(&status,MPI_INT,&count);
						recvindex[index].resize(count);
						MPI_Recv(&recvindex[index][0],count,MPI_INT,procorder[i],0,
							MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
					if(bmerit[index].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
						MPI_Send(&bmerit[index][0],(int)bmerit[index].size(),
						MPI_DOUBLE,procorder[i],0,MPI_COMM_WORLD);
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					if(status.MPI_TAG==1)
						MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					else
					{
						int count;
						MPI_Get_count(&status,MPI_DOUBLE,&count);
						recvmerit[index].resize(count);
						MPI_Recv(&recvmerit[index][0],count,MPI_DOUBLE,procorder[i],0,
							MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
				}
				else
				{
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					if(status.MPI_TAG==1)
						MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					else
					{
						int count;
						MPI_Get_count(&status,MPI_INT,&count);
						recvindex[index].resize(count);
						MPI_Recv(&recvindex[index][0],count,MPI_INT,procorder[i],0,
							MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
					if(bremove[index].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
						MPI_Send(&bremove[index][0],(int)bremove[index].size(),
						MPI_INT,procorder[i],0,MPI_COMM_WORLD);
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					if(status.MPI_TAG==1)
						MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					else
					{
						int count;
						MPI_Get_count(&status,MPI_DOUBLE,&count);
						recvmerit[index].resize(count);
						MPI_Recv(&recvmerit[index][0],count,MPI_DOUBLE,procorder[i],0,
							MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
					if(bmerit[index].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
						MPI_Send(&bmerit[index][0],(int)bmerit[index].size(),
						MPI_DOUBLE,procorder[i],0,MPI_COMM_WORLD);
				}
			}
		}

		vector<int> res,bad;
		vector<vector<int> > const& Nghost=tess.GetGhostIndeces();
		vector<vector<int> > const& DupPoints=tess.GetDuplicatedPoints();
		for(int i=0;i<nlist;++i)
		{
			if(recvindex[i].empty()||sentpoints[i].empty())
				continue;
			vector<int> Nghostindex(recvindex[i].size());
			for(int j=0;j<(int)Nghostindex.size();++j)
				Nghostindex[j]=Nghost[i][recvindex[i][j]];
			vector<int> indeces;
			sort_index(Nghostindex,indeces);
			sort(Nghostindex.begin(),Nghostindex.end());
			for(int j=0;j<(int)bremove[i].size();++j)
			{
				vector<int> neigh=tess.GetNeighbors(DupPoints[i][bremove[i][j]]);
				for(int k=0;k<(int)neigh.size();++k)
				{
					if(!binary_search(Nghostindex.begin(),Nghostindex.end(),neigh[k]))
						continue;
					int index=lower_bound(Nghostindex.begin(),Nghostindex.end(),neigh[k])-
						Nghostindex.begin();
					if(index<(int)Nghostindex.size())
						if(recvmerit[i][indeces[index]]>=bmerit[i][j])
							bad.push_back(DupPoints[i][bremove[i][j]]);
				}
			}
		}
		sort(bad.begin(),bad.end());
		bad=unique(bad);
		res=RemoveList(toremove,bad);
		return res;
}

#endif

void PeriodicUpdateCells(vector<Primitive> &cells,vector<vector<double> > &tracers,
	vector<size_t> &customevolutions,vector<vector<int> > const& sentcells,
	int npoints)
{
	if(sentcells.empty())
		return;
	bool traceractive=!tracers.empty();
	if(traceractive)
		traceractive=!tracers[0].empty();
	int n=(int)sentcells.size();
	int totalsent=0;
	for(int i=0;i<n;++i)
		if(!sentcells[i].empty())
			totalsent+=(int)sentcells[i].size();
	cells.resize(npoints-totalsent);
	customevolutions.resize(npoints-totalsent);
	if(traceractive)
		tracers.resize(npoints-totalsent);
	for(int i=0;i<n;++i)
	{
		if(sentcells[i].empty())
			continue;
		vector<Primitive> ptemp=VectorValues(cells,sentcells[i]);
		vector<size_t> ctemp=VectorValues(customevolutions,sentcells[i]);
		if(!ptemp.empty())
			cells.insert(cells.end(),ptemp.begin(),ptemp.end());
		if(!ctemp.empty())
			customevolutions.insert(customevolutions.end(),ctemp.begin(),ctemp.end());
		if(traceractive)
		{
			vector<vector<double> > temp=VectorValues(tracers,sentcells[i]);
			if(!temp.empty())
				tracers.insert(tracers.end(),temp.begin(),temp.end());
		}
	}
}

void PeriodicVelocityExchange(vector<Vector2D> &vel,
	vector<vector<int> > const& sentcells,int npoints)
{
	int n=(int)sentcells.size();
	int totalsent=0;
	for(int i=0;i<n;++i)
		if(!sentcells[i].empty())
			totalsent+=(int)sentcells[i].size();
	vel.resize(npoints-totalsent);

	for(int i=0;i<n;++i)
	{
		if(sentcells[i].empty())
			continue;
		vector<Vector2D> temp=VectorValues(vel,sentcells[i]);
		if(!temp.empty())
			vel.insert(vel.end(),temp.begin(),temp.end());
	}
}

void PeriodicGradExchange(vector<ReducedPrimitiveGradient2D> &grad,
	vector<vector<int> > const& sentcells,int npoints)
{
	int n=(int)sentcells.size();
	int totalsent=0;
	for(int i=0;i<n;++i)
		if(!sentcells[i].empty())
			totalsent+=(int)sentcells[i].size();
	grad.resize(npoints-totalsent);

	for(int i=0;i<n;++i)
	{
		if(sentcells[i].empty())
			continue;
		vector<ReducedPrimitiveGradient2D> ptemp=VectorValues(grad,sentcells[i]);
		if(!ptemp.empty())
			grad.insert(grad.end(),ptemp.begin(),ptemp.end());
	}
}
