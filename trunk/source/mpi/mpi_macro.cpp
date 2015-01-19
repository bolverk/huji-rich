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
		int nlist=static_cast<int>(sentprocs.size());
		const int rank = get_mpi_rank();
		const int worldsize = get_mpi_size();
		const vector<int> procorder=GetProcOrder(rank,worldsize);
		int n=worldsize-1;
		padd.resize(static_cast<size_t>(nlist));
		tadd.resize(static_cast<size_t>(nlist));
		// Send the data
		for(size_t i=0;i<static_cast<size_t>(n);++i)
		{
		  int index=static_cast<int>(Find(sentprocs.begin(),sentprocs.end(),procorder[i])
					     -sentprocs.begin());
			if(index<nlist)
			{
			  if(rank<procorder[i])
				{
				  vector<Primitive> ptemp=VectorValues(cells,sentcells[static_cast<size_t>(index)]);
				  MPI_SendVectorPrimitive(ptemp,procorder[i],0,MPI_COMM_WORLD);
				  MPI_RecvVectorPrimitive(padd[static_cast<size_t>(index)],procorder[i],0,MPI_COMM_WORLD,eos);
					if(traceractive)
					{
					  vector<vector<double> > ttemp=VectorValues(tracers,sentcells[static_cast<size_t>(index)]);
					  MPI_SendVectorTracer(ttemp,procorder[i],0,MPI_COMM_WORLD);
					  MPI_RecvVectorTracer(tadd[static_cast<size_t>(index)],procorder[i],0,MPI_COMM_WORLD,
							       static_cast<int>(tracers[0].size()));
					}
				}
				else
				{
				  MPI_RecvVectorPrimitive(padd[static_cast<size_t>(index)],procorder[i],0,MPI_COMM_WORLD,eos);
				  vector<Primitive> ptemp=VectorValues(cells,sentcells[static_cast<size_t>(index)]);
				  MPI_SendVectorPrimitive(ptemp,procorder[i],0,MPI_COMM_WORLD);
					if(traceractive)
					{
					  MPI_RecvVectorTracer(tadd[static_cast<size_t>(index)],procorder[i],0,MPI_COMM_WORLD,
							       static_cast<int>(tracers[0].size()));
					  vector<vector<double> > ttemp=VectorValues(tracers,sentcells[static_cast<size_t>(index)]);
					  MPI_SendVectorTracer(ttemp,procorder[i],0,MPI_COMM_WORLD);
					}
				}
			}
		}
	}

	void GradVectorToDouble(vector<ReducedPrimitiveGradient2D> const& vec,
		vector<double> &res)
	{
		int n=static_cast<int>(vec.size());
		int gradlength=8;
		if(!vec[0].tracers.empty())
		  gradlength+=static_cast<int>(vec[0].tracers.size());
		if(n*gradlength!=static_cast<int>(res.size()))
		{
			UniversalError eo("Sizes do not match in GradVectorToDouble");
			throw eo;
		}
		for(size_t i=0;i<static_cast<size_t>(n);++i)
		{
			res[static_cast<size_t>(gradlength)*i]=vec[i].density.x;
			res[static_cast<size_t>(gradlength)*i+1]=vec[i].density.y;
			res[static_cast<size_t>(gradlength)*i+2]=vec[i].pressure.x;
			res[static_cast<size_t>(gradlength)*i+3]=vec[i].pressure.y;
			res[static_cast<size_t>(gradlength)*i+4]=vec[i].xvelocity.x;
			res[static_cast<size_t>(gradlength)*i+5]=vec[i].xvelocity.y;
			res[static_cast<size_t>(gradlength)*i+6]=vec[i].yvelocity.x;
			res[static_cast<size_t>(gradlength)*i+7]=vec[i].yvelocity.y;
			for(size_t j=0;j<static_cast<size_t>(gradlength-8)/2;++j)
			{
				res[static_cast<size_t>(gradlength)*i+j*2+8]=vec[i].tracers[j].x;
				res[static_cast<size_t>(gradlength)*i+j*2+9]=vec[i].tracers[j].y;
			}
		}
	}

	void DoubleVectorToGrad(vector<ReducedPrimitiveGradient2D> &vec,
		vector<double> const& temp,int gradlength)
	{
	  int n=static_cast<int>(temp.size())/gradlength;
		if(n!=static_cast<int>(vec.size()))
		{
			UniversalError eo("Sizes do not match in DoubleVectorToGrad");
			throw eo;
		}
		for(size_t i=0;i<static_cast<size_t>(n);++i)
		{
			vec[i].density.x=temp[static_cast<size_t>(gradlength)*i];
			vec[i].density.y=temp[static_cast<size_t>(gradlength)*i+1];
			vec[i].pressure.x=temp[static_cast<size_t>(gradlength)*i+2];
			vec[i].pressure.y=temp[static_cast<size_t>(gradlength)*i+3];
			vec[i].xvelocity.x=temp[static_cast<size_t>(gradlength)*i+4];
			vec[i].xvelocity.y=temp[static_cast<size_t>(gradlength)*i+5];
			vec[i].yvelocity.x=temp[static_cast<size_t>(gradlength)*i+6];
			vec[i].yvelocity.y=temp[static_cast<size_t>(gradlength)*i+7];
			for(size_t j=0;j<static_cast<size_t>(gradlength-8)/2;++j)
			{
				Vector2D vtemp(temp[static_cast<size_t>(gradlength)*i+2*j+8],temp[static_cast<size_t>(gradlength)*i+2*j+9]);
				vec[i].tracers.push_back(vtemp);
			}
		}
	}

	void DoubleVectorToTracer(vector<vector<double> > &tracer,vector<double>
		const& data,int tracerlength)
	{
	  int n=static_cast<int>(data.size())/tracerlength;
		if(n!=static_cast<int>(tracer.size()))
		{
			UniversalError eo("Sizes do not match in DoubleVectorToTracer");
			throw eo;
		}
		for(size_t i=0;i<static_cast<size_t>(n);++i)
		{
		  tracer[i].resize(static_cast<size_t>(tracerlength));
		  for(size_t j=0;j<static_cast<size_t>(tracerlength);++j)
		    tracer[i][j]=data[i*static_cast<size_t>(tracerlength)+j];
		}
	}

	void TracerVectorToDouble(vector<vector<double> > const& ttemp,vector<double>
		&dtracer)
	{
	  int tracerlength=static_cast<int>(ttemp[0].size());
		int n=static_cast<int>(ttemp.size());
		if(n*tracerlength!=static_cast<int>(dtracer.size()))
		{
			UniversalError eo("Sizes do not match in TracerVectorToDouble");
			throw eo;
		}
		for(size_t i=0;i<static_cast<size_t>(n);++i)
		  for(size_t j=0;j<static_cast<size_t>(tracerlength);++j)
		    dtracer[i*static_cast<size_t>(tracerlength)+j]=ttemp[i][j];
	}

	void DoubleVectorToPrimitve(vector<double> const& temp,vector<Primitive> &vec,
		EquationOfState const& eos)
	{
		int ntotal=static_cast<int>(temp.size()/4);
		if(ntotal!=static_cast<int>(vec.size()))
		{
			UniversalError eo("Sizes do not match in DoubleVectorToPrimitve");
			throw eo;
		}
		for(size_t i=0;i<static_cast<size_t>(ntotal);++i)
		{
			Primitive ptemp(temp[4*i],temp[4*i+1],Vector2D(temp[4*i+2],temp[4*i+3]),0,0);
			ptemp.Energy=eos.dp2e(ptemp.Density,ptemp.Pressure);
			ptemp.SoundSpeed=eos.dp2c(ptemp.Density,ptemp.Pressure);
			vec[i]=ptemp;
		}
	}

	void PrimitiveVectorToDouble(vector<Primitive> const& vec,vector<double> &res)
	{
		int n=static_cast<int>(vec.size());
		if(n*4!=static_cast<int>(res.size()))
		{
			UniversalError eo("Sizes do not match in PrimitiveVectorToDouble");
			throw eo;
		}
		for(size_t i=0;i<static_cast<size_t>(n);++i)
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
	  int n=static_cast<int>(data.size()/2);
		if(n!=static_cast<int>(vec.size()))
		{
			UniversalError eo("Sizes do not match in DoubleVectorToVector2D");
			throw eo;
		}
		for(size_t i=0;i<static_cast<size_t>(n);++i)
			vec[i]=Vector2D(data[2*i],data[2*i+1]);
	}

	void Vector2DVectorToDouble(vector<double> &res,vector<Vector2D> const& data)
	{
		int n=static_cast<int>(data.size());
		if(n*2!=static_cast<int>(res.size()))
		{
			UniversalError eo("Sizes do not match in Vector2DVectorToDouble");
			throw eo;
		}
		for(size_t i=0;i<static_cast<size_t>(n);++i)
		{
			res[2*i]=data[i].x;
			res[2*i+1]=data[i].y;
		}
	}

	int FindLoc(vector<int> const& vec,int data,int occur)
	{
		for(size_t i=0;i<vec.size();++i)
		{
		  if(vec[i]==data){
				if(occur==0)
				  return static_cast<int>(i);
				else
					--occur;
		  }
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
	int n=static_cast<int>(vec.size());
	if(vec.empty())
	{
		double temp=0;
		int err=MPI_Send(&temp,1,MPI_DOUBLE,dest,1,comm);
		return err;
	}
	vector<double> temp(static_cast<size_t>(n)*2);
	for(size_t i=0;i<static_cast<size_t>(n);++i)
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
	vector<double> temp(static_cast<size_t>(n));
	vec.resize(static_cast<size_t>(n)/2);
	n/=2;
	int err=0;
	err=MPI_Recv(&temp[0],2*n,MPI_DOUBLE,source,tag,comm,MPI_STATUS_IGNORE);
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
		vec[i].x=temp[2*i];
		vec[i].y=temp[2*i+1];
	}
	return err;
}

int MPI_VectorBcast_Vector2D(vector<Vector2D> &vec,int root, MPI_Comm comm,int rank)
{
	int n=static_cast<int>(vec.size());
	vector<double> temp(static_cast<size_t>(n)*2);
	if(rank==root)
	{
	  for(size_t i=0;i<static_cast<size_t>(n);++i)
		{
			temp[2*i]=vec[i].x;
			temp[2*i+1]=vec[i].y;
		}
	}
	int err=MPI_Bcast(&temp[0],n*2,MPI_DOUBLE,root,comm);
	if(rank!=root)
	{
	  for(size_t i=0;i<static_cast<size_t>(n);++i)
		{
			vec[i].x=temp[2*i];
			vec[i].y=temp[2*i+1];
		}
	}
	return err;
}

bool PointInsideCell(Tessellation const& tess,int cell_index,Vector2D const & point)
{
	vector<Vector2D> cpoints;
	ConvexHull(cpoints, &tess, cell_index);
	boost::array<Vector2D, 3> tocheck;
	tocheck[2] = point;
	for (size_t i = 0; i<cpoints.size(); ++i)
	{
		tocheck[0] = cpoints[i];
		tocheck[1] = cpoints[(i + 1) % cpoints.size()];
		if (orient2d(TripleConstRef<Vector2D>
			     (cpoints[i],
			      cpoints[(i+1)%cpoints.size()],
			      point))<0)
			return false;
	}
	return true;
}

void ConvertDoubleToVector2D(vector<Vector2D> & res,vector<double> const& vec)
{
	res.resize(vec.size()/2);
	for(size_t i=0;i<res.size();++i)
		res[i].Set(vec[2*i],vec[2*i+1]);
}

void ConvertVector2DToDouble(vector<Vector2D> const& vec,vector<double> &res)
{
  int n=static_cast<int>(vec.size());
  res.resize(static_cast<size_t>(n)*2);
  for(size_t i=0;i<static_cast<size_t>(n);++i)
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
	int n=static_cast<int>(procorder.size());
	int nlist=static_cast<int>(proclist.size());
	const int rank = get_mpi_rank();
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  int index=static_cast<int>(find(proclist.begin(),proclist.end(),procorder[i])-proclist.begin());
		// Do we talk with this processor?
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
			  if(tosend[static_cast<size_t>(index)].empty())
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD);
				}
				else
				  MPI_VectorSend_Vector2D(tosend[static_cast<size_t>(index)],procorder[i],0,
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
				if(tosend[static_cast<size_t>(index)].empty())
				{
					double temp=0;
					MPI_Send(&temp,1,MPI_DOUBLE,procorder[i],1,MPI_COMM_WORLD);
				}
				else
				  MPI_VectorSend_Vector2D(tosend[static_cast<size_t>(index)],procorder[i],0,
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
	vector<int> procorder(static_cast<size_t>(worldsize)-1);
	if(rank==0)
	{
		for(size_t i=0;i<static_cast<size_t>(worldsize-1);++i)
		  procorder[i]=static_cast<int>(i)+1;
	}
	if(rank!=(worldsize-1))
	{
		for(size_t i=0;i<static_cast<size_t>(worldsize-1);++i)
		{
		  int temp=(static_cast<int>(i)-rank+worldsize)%(worldsize-1);
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
			int half=static_cast<int>(worldsize/2);
			for(size_t i=0;i<static_cast<size_t>(worldsize-1);++i)
			{
				if(i%2==0)
				  procorder[i]=half+static_cast<int>(i)/2;
				else
					procorder[i]=static_cast<int>(i/2)+1;
			}
			procorder[static_cast<size_t>(worldsize)-2]=0;
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
		cei_to_send_(cei_to_send),ntracer_(ntracer), 
		c_received_(),
		t_received_(), 
		cei_received_() {}

		void sendInfo(int address)
		{
			MPI_SendVectorConserved(c_to_send_,address,0,MPI_COMM_WORLD);
			MPI_SendVectorTracer(t_to_send_,address,0,MPI_COMM_WORLD);
			if(!cei_to_send_.empty())
			{
				vector<unsigned> buf = mass_static_cast<unsigned,size_t>(cei_to_send_);
				MPI_Send(&buf[0],static_cast<int>(buf.size()),MPI_UNSIGNED,address,0,MPI_COMM_WORLD);
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
				vector<unsigned> buf(static_cast<size_t>(count));
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
	  ntracer=static_cast<int>(tracers[0].size());

	// Take care of self send hydro
	int nlist=static_cast<int>(sentprocs.size());
	ptoadd.clear();
	ttoadd.clear();
	ctoadd.clear();

	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);

	int n=worldsize-1;
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  int index=static_cast<int>(find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin());
		// Do we talk with this processor?
		if(index<nlist)
		{
		  ExtensiveCommunicator extensive_communicator(VectorValues(cons,sentcells[static_cast<size_t>(index)]),
							       VectorValues(tracers,sentcells[static_cast<size_t>(index)]),
							       VectorValues(customevolutions,sentcells[static_cast<size_t>(index)]),ntracer);
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
	int nlist=static_cast<int>(sentprocs.size());
	btoadd.clear();

	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);

	int n=worldsize-1;
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  int index=static_cast<int>(find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin());
		// Do we talk with this processor?
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<char> btemp=VectorValues(shockedcells,sentcells[static_cast<size_t>(index)]);
				if(!btemp.empty())
					MPI_Send(&btemp[0],static_cast<int>(btemp.size()),MPI_CHAR,procorder[i],0,
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
					vector<char> toadd(static_cast<size_t>(count));
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
					vector<char> toadd(static_cast<size_t>(count));
					MPI_Recv(&toadd[0],count,MPI_CHAR,procorder[i],0,
						MPI_COMM_WORLD,&status);
					btoadd.insert(btoadd.end(),toadd.begin(),toadd.end());
				}
				vector<char> btemp=VectorValues(shockedcells,sentcells[static_cast<size_t>(index)]);
				if(!btemp.empty())
					MPI_Send(&btemp[0],static_cast<int>(btemp.size()),MPI_CHAR,procorder[i],0,
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
	int nlist=static_cast<int>(sentprocs.size());
	toadd.clear();

	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);

	int n=worldsize-1;
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  int index=static_cast<int>(find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin());
		// Do we talk with this processor?
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<double> temp=VectorValues(vec,sentcells[static_cast<size_t>(index)]);
				if(!temp.empty())
					MPI_Send(&temp[0],static_cast<int>(temp.size()),MPI_DOUBLE,procorder[i],0,
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
					vector<double> ttoadd(static_cast<size_t>(count));
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
					vector<double> ttoadd(static_cast<size_t>(count));
					MPI_Recv(&ttoadd[0],count,MPI_DOUBLE,procorder[i],0,
						MPI_COMM_WORLD,&status);
					toadd.insert(toadd.end(),ttoadd.begin(),ttoadd.end());
				}
				vector<double> temp=VectorValues(vec,sentcells[static_cast<size_t>(index)]);
				if(!temp.empty())
					MPI_Send(&temp[0],static_cast<int>(temp.size()),MPI_DOUBLE,procorder[i],0,
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
	int nlist=static_cast<int>(sentprocs.size());
	toadd.clear();

	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);

	int n=worldsize-1;
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  int index=static_cast<int>(find(sentprocs.begin(),sentprocs.end(),procorder[i])-sentprocs.begin());
		// Do we talk with this processor?
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<Vector2D> vtemp=VectorValues(points,sentcells[static_cast<size_t>(index)]);
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
				vtemp=VectorValues(points,sentcells[static_cast<size_t>(index)]);
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
			  vector<unsigned> buf(static_cast<size_t>(count));
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
	int nlist=static_cast<int>(sentprocs.size());
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int n=worldsize-1;
	vector<vector<Primitive> > padd(static_cast<size_t>(nlist));
	vector<vector<size_t> > cadd(static_cast<size_t>(nlist));

	// Send the data
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  int index=static_cast<int>(Find(sentprocs.begin(),sentprocs.end(),procorder[i])
				     -sentprocs.begin());
		if(index==nlist)
			continue;
		HydroCommunicator hydro_communicator(eos,
			VectorValues(cells,sentcells[static_cast<size_t>(index)]),
			VectorValues(customevolutions,sentcells[static_cast<size_t>(index)]));
		marshal_communication(hydro_communicator,procorder[i],rank<procorder[i]);
		padd[static_cast<size_t>(index)] = hydro_communicator.getReceivedPrimitives();
		cadd[static_cast<size_t>(index)] = hydro_communicator.getReceivedCustomEvolutionIndices();
	}
	// ReArrange the data
	cells.resize(static_cast<size_t>(totalpoints));
	customevolutions.resize(static_cast<size_t>(totalpoints));
	for(size_t i=0;i<static_cast<size_t>(nlist);++i)
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
	int nlist=static_cast<int>(sentprocs.size());
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	vector<vector<vector<double> > > tadd(static_cast<size_t>(nlist));
	// Send the data
	for(size_t i=0;i<static_cast<size_t>(worldsize)-1;++i)
	{
	  int index=static_cast<int>(Find(sentprocs.begin(),sentprocs.end(),procorder[i])
				     -sentprocs.begin());
		if(index<nlist)
		{
			TracerCommunicator tc(VectorValues(tracers,sentcells[static_cast<size_t>(index)]),
					      static_cast<int>(tracers[0].size()));
			marshal_communication(tc,procorder[i],
				rank<procorder[i]);
			tadd[static_cast<size_t>(index)] = tc.getReply();
		}
	}
	// ReArrange the data
	tracers.resize(static_cast<size_t>(totalpoints));
	for(size_t i=0;i<static_cast<size_t>(nlist);++i)
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
	int n=static_cast<int>(vec.size());
	vector<double> tosend(n*4);
	for(size_t i=0;i<static_cast<size_t>(n);++i)
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
	vector<double> temp(static_cast<size_t>(nrecv));
	MPI_Recv(&temp[0],nrecv,MPI_DOUBLE,dest,tag,comm,&status);
	int ntotal=nrecv/4;
	vec.reserve(ntotal);
	for(size_t i=0;i<static_cast<size_t>(ntotal);++i)
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
	int n=static_cast<int>(vec.size());
	int ntracer=static_cast<int>(vec[0].size());
	vector<double> tosend(n*ntracer);
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  for(size_t j=0;j<static_cast<size_t>(ntracer);++j)
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
	vector<double> temp(static_cast<size_t>(nrecv));
	MPI_Recv(&temp[0],nrecv,MPI_DOUBLE,dest,MPI_ANY_TAG,comm,&status);
	int length=nrecv/ntracer;
	vec.resize(static_cast<size_t>(length));
	for(size_t i=0;i<static_cast<size_t>(length);++i)
	{
	  vec[i].resize(static_cast<size_t>(ntracer));
		for(size_t j=0;j<static_cast<size_t>(ntracer);++j)
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
	int n=static_cast<int>(vec.size());
	int gradlength=2*static_cast<int>(vec[0].tracers.size())+8;
	vector<double> tosend(n*gradlength);
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
		tosend[static_cast<size_t>(gradlength)*i]=vec[i].density.x;
		tosend[static_cast<size_t>(gradlength)*i+1]=vec[i].density.y;
		tosend[static_cast<size_t>(gradlength)*i+2]=vec[i].pressure.x;
		tosend[static_cast<size_t>(gradlength)*i+3]=vec[i].pressure.y;
		tosend[static_cast<size_t>(gradlength)*i+4]=vec[i].xvelocity.x;
		tosend[static_cast<size_t>(gradlength)*i+5]=vec[i].xvelocity.y;
		tosend[static_cast<size_t>(gradlength)*i+6]=vec[i].yvelocity.x;
		tosend[static_cast<size_t>(gradlength)*i+7]=vec[i].yvelocity.y;
		for(size_t j=0;j<static_cast<size_t>(gradlength-8)/2;++j)
		{
			tosend[static_cast<size_t>(gradlength)*i+j*2+8]=vec[i].tracers[j].x;
			tosend[static_cast<size_t>(gradlength)*i+j*2+9]=vec[i].tracers[j].y;
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
	vector<double> temp(static_cast<size_t>(ntotal));
	MPI_Recv(&temp[0],ntotal,MPI_DOUBLE,dest,tag,comm,&status);
	int n=ntotal/gradlength;
	vec.resize(static_cast<size_t>(n));
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
		vec[i].density.x=temp[static_cast<size_t>(gradlength)*i];
		vec[i].density.y=temp[static_cast<size_t>(gradlength)*i+1];
		vec[i].pressure.x=temp[static_cast<size_t>(gradlength)*i+2];
		vec[i].pressure.y=temp[static_cast<size_t>(gradlength)*i+3];
		vec[i].xvelocity.x=temp[static_cast<size_t>(gradlength)*i+4];
		vec[i].xvelocity.y=temp[static_cast<size_t>(gradlength)*i+5];
		vec[i].yvelocity.x=temp[static_cast<size_t>(gradlength)*i+6];
		vec[i].yvelocity.y=temp[static_cast<size_t>(gradlength)*i+7];
		for(size_t j=0;j<static_cast<size_t>(gradlength-8)/2;++j)
		{
			Vector2D vtemp(temp[static_cast<size_t>(gradlength)*i+2*j+8],temp[static_cast<size_t>(gradlength)*i+2*j+9]);
			vec[i].tracers.push_back(vtemp);
		}
	}
}

void SendRecvVelocity(vector<Vector2D> &vel,vector<vector<int> >const& sentcells,
	vector<int> sentprocs,vector<vector<int> > const& Nghost,int totalpoints)
{
	int nlist=static_cast<int>(sentprocs.size());
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int n=worldsize-1;
	vector<vector<Vector2D> > tadd(static_cast<size_t>(nlist));
	// Send the data
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  int index=static_cast<int>(Find(sentprocs.begin(),sentprocs.end(),procorder[i])
				     -sentprocs.begin());
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<Vector2D> ttemp=VectorValues(vel,sentcells[static_cast<size_t>(index)]);
				MPI_VectorSend_Vector2D(ttemp,procorder[i],0,MPI_COMM_WORLD);
				MPI_VectorRecv_Vector2D(tadd[static_cast<size_t>(index)],procorder[i],0,MPI_COMM_WORLD);
			}
			else
			{
				MPI_VectorRecv_Vector2D(tadd[static_cast<size_t>(index)],procorder[i],0,MPI_COMM_WORLD);
				vector<Vector2D> ttemp=VectorValues(vel,sentcells[static_cast<size_t>(index)]);
				MPI_VectorSend_Vector2D(ttemp,procorder[i],0,MPI_COMM_WORLD);
			}
		}
	}
	// ReArrange the data
	vel.resize(static_cast<size_t>(totalpoints));
	for(size_t i=0;i<static_cast<size_t>(nlist);++i)
		ListExchange(vel,Nghost[i],tadd[i]);
}

void SendRecvGrad(vector<ReducedPrimitiveGradient2D> &grads,
	vector<vector<int> >const& sentcells,vector<int> sentprocs,
	vector<vector<int> > const& Nghost,int totalpoints)
{
	if(grads.empty())
		return;
	int nlist=static_cast<int>(sentprocs.size());
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int n=worldsize-1;
	vector<vector<ReducedPrimitiveGradient2D> > tadd(static_cast<size_t>(nlist));
	int gradlength=8+2*static_cast<int>(grads[0].tracers.size());
	// Send the data
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  int index=static_cast<int>(Find(sentprocs.begin(),sentprocs.end(),procorder[i])
				     -sentprocs.begin());
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				vector<ReducedPrimitiveGradient2D> ttemp=VectorValues(grads,sentcells[static_cast<size_t>(index)]);
				MPI_SendVectorGrad(ttemp,procorder[i],0,MPI_COMM_WORLD);
				MPI_RecvVectorGrad(tadd[static_cast<size_t>(index)],procorder[i],0,MPI_COMM_WORLD,gradlength);
			}
			else
			{
				vector<ReducedPrimitiveGradient2D> ttemp=VectorValues(grads,sentcells[static_cast<size_t>(index)]);
				MPI_RecvVectorGrad(tadd[static_cast<size_t>(index)],procorder[i],0,MPI_COMM_WORLD,gradlength);
				MPI_SendVectorGrad(ttemp,procorder[i],0,MPI_COMM_WORLD);
			}
		}
	}
	// ReArrange the data
	grads.resize(static_cast<size_t>(totalpoints));
	for(size_t i=0;i<static_cast<size_t>(nlist);++i)
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
	int n=static_cast<int>(vec.size());
	vector<double> tosend(n*4);
	for(size_t i=0;i<static_cast<size_t>(n);++i)
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
	vector<double> temp(static_cast<size_t>(nrecv));
	MPI_Recv(&temp[0],nrecv,MPI_DOUBLE,dest,tag,comm,&status);
	int ntotal=nrecv/4;
	vec.reserve(ntotal);
	for(size_t i=0;i<static_cast<size_t>(ntotal);++i)
	{
		Conserved ctemp(temp[4*i],Vector2D(temp[4*i+2],temp[4*i+3]),temp[4*i+1]);
		vec.push_back(ctemp);
	}
}

void SendRecvGhostIndeces(vector<vector<int> > &GhostIndeces,vector<int>
	const& BoundaryPoints,vector<vector<int> > const& SentPoints,vector<int> const&
	SentProcs)
{
	int nprocs=static_cast<int>(SentProcs.size());
	int nbound=static_cast<int>(BoundaryPoints.size());
	/*
	const int rank = get_mpi_rank();
	const int ws = get_mpi_size();
	*/
	vector<MPI_Status> status(static_cast<size_t>(nprocs));
	vector<MPI_Request> req(static_cast<size_t>(nprocs));

	vector<vector<int> > tosend(static_cast<size_t>(nprocs)),torecv(static_cast<size_t>(nprocs));
	vector<int> sentme,flags(nprocs,0);
	int itemp=-1;

	for(size_t i=0;i<static_cast<size_t>(nprocs);++i)
	{
		// Send the data
		if(!SentPoints[i].empty())
		{
			// Find the relevant points
			vector<int> temp(SentPoints[i]);
			sort(temp.begin(),temp.end());
			for(size_t j=0;j<static_cast<size_t>(nbound);++j)
			{
			  const int index=static_cast<int>(lower_bound(temp.begin(),temp.end(),BoundaryPoints[j])
							   -temp.begin());
				if(index<static_cast<int>(temp.size()))
					tosend[i].push_back(index);
			}
			MPI_Isend(&tosend[i][0],static_cast<int>(tosend.size()),MPI_INT,SentProcs[i],0,
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
	  for(size_t i=0;i<static_cast<size_t>(nprocs);++i)
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
						torecv[i].resize(static_cast<size_t>(count));
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
	GhostIndeces.resize(static_cast<size_t>(nprocs));
	for(size_t i=0;i<static_cast<size_t>(nprocs);++i)
	{
	  int loc=static_cast<int>(lower_bound(occur.begin(),occur.end(),sentme[i])-occur.begin());
		int index=FindLoc(SentProcs,sentme[i],occur[loc]);
		++occur[loc];
		GhostIndeces[static_cast<size_t>(index)]=torecv[i];
	}
}

namespace
{
	// Nghost is the NghostIndex
	// ghost is the index in cor of the ghost removed points
	// Sentindex is the index of the proc we are dealing with
	// returns all of the local neighboring points of the removed ghost point
  vector<vector<int> > FindLocalNeighbors(vector<int> const& /*Nghost*/,vector<int>  &ghost,
		int SentIndex,Tessellation const& tess)
	{
		vector<vector<int> > const& duplicated=tess.GetDuplicatedPoints();
		vector<int> const& Sent=duplicated[SentIndex];
		vector<int> index;
		sort_index(ghost,index);
		sort(ghost.begin(),ghost.end());
		int nreal=static_cast<int>(ghost.size());
		vector<vector<int> > res(ghost.size());
		for(size_t i=0;i<Sent.size();++i)
		{ //do i need to check from other procs??????
			vector<int> neigh=tess.GetNeighbors(Sent[i]);
			for(size_t j=0;j<neigh.size();++j)
			{
				if(binary_search(ghost.begin(),ghost.end(),neigh[j]))
				{
				  const size_t index2=static_cast<size_t>(lower_bound(ghost.begin(),ghost.end(),neigh[j])-
									  ghost.begin());
				  res[index[index2]].push_back(Sent[i]);
				}
			}
		}
		for(size_t i=0;i<static_cast<size_t>(nreal);++i)
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
		int nprocs=static_cast<int>(BoundaryRemove.size());
		vector<vector<int> > const& Nghost=tess.GetGhostIndeces();
		localneigh.clear();
		ghostneigh.clear();
		for(size_t i=0;i<static_cast<size_t>(nprocs);++i)
		{
			if(BoundaryRemove.empty())
				continue;
			vector<int> ghosts(BoundaryRemove[i].size());
			for(size_t j=0;j<BoundaryRemove[i].size();++j)
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
			vector<vector<int> > temp=FindLocalNeighbors(Nghost[i],ghosts,static_cast<int>(i),tess);
			for(size_t j=0;j<BoundaryRemove[i].size();++j)
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
	int nlist=static_cast<int>(proclist.size());
	vector<vector<int> > recv(static_cast<size_t>(nlist));
	// Talk with other procs
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int n=worldsize-1;
	int temp;
	// Send the data
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
	  int index=static_cast<int>(Find(proclist.begin(),proclist.end(),procorder[i])
				     -proclist.begin());
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				if(ToSend[static_cast<size_t>(index)].empty())
					MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
				  MPI_Send(&ToSend[static_cast<size_t>(index)][0],static_cast<int>(ToSend[static_cast<size_t>(index)].size()),MPI_INT,
					procorder[i],0,MPI_COMM_WORLD);
				MPI_Status stat;
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					recv[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
					MPI_Recv(&recv[static_cast<size_t>(index)][0],count,MPI_INT,procorder[i],0,
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
					recv[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
					MPI_Recv(&recv[static_cast<size_t>(index)][0],count,MPI_INT,procorder[i],0,
						MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if(ToSend[static_cast<size_t>(index)].empty())
					MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
				  MPI_Send(&ToSend[static_cast<size_t>(index)][0],static_cast<int>(ToSend[static_cast<size_t>(index)].size()),MPI_INT,
					procorder[i],0,MPI_COMM_WORLD);
			}
		}
	}
	vector<vector<int> > ToSend2(static_cast<size_t>(nlist));
	for(size_t i=0;i<static_cast<size_t>(nlist);++i)
	  for(size_t j=0;j<recv[i].size();++j)
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
	for(size_t i=0;i<static_cast<size_t>(nlist);++i)
	{
		if(!ToSend[i].empty())
		{
		  for(size_t j=0;j<ToSend[i].size();++j)
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
	for(int k=0;k<static_cast<int>(ToRemove.size());++k)
	{
	  for(size_t i=0;i<static_cast<size_t>(nlist);++i)
		{
			if(binary_search(ToSend[i].begin(),ToSend[i].end(),ToRemove[k]))
			{
			  int index2=static_cast<int>(lower_bound(ToSend[i].begin(),ToSend[i].end(),
								  ToRemove[k])-ToSend[i].begin());
				rescells.push_back(padd[i][index2]);
				if(traceractive)
					restracer.push_back(tadd[i][index2]);
				break;
			}
		}
		if(static_cast<int>(rescells.size())!=(k+1))
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
	int nlist=static_cast<int>(sentprocs.size());
	int n=worldsize-1;
	vector<int> proclist=tess.GetDuplicatedProcs();
	vector<vector<int> > RecvPoints(proclist.size());
	vector<vector<vector<int> > > RecvNeigh(proclist.size());
	int trash,trashr;
	MPI_Status stat;
	for(size_t i=0;i<static_cast<size_t>(n);++i)
	{
		//	cout<<"rank "<<rank<<" i "<<i<<" talk with "<<procorder[i]<<endl;
	  int index=static_cast<int>(Find(sentprocs.begin(),sentprocs.end(),procorder[i])
				     -sentprocs.begin());
		if(index<nlist)
		{
			if(rank<procorder[i])
			{
				if(BoundaryRemove[static_cast<size_t>(index)].empty())
					MPI_Send(&trash,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
				  MPI_Send(&BoundaryRemove[static_cast<size_t>(index)][0],static_cast<int>(BoundaryRemove[static_cast<size_t>(index)].size()),
					MPI_INT,procorder[i],0,MPI_COMM_WORLD);
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&trashr,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					RecvPoints[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
					MPI_Recv(&RecvPoints[static_cast<size_t>(index)][0],count,MPI_INT,procorder[i],0,
						MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				// create the lengths of neighbors
				vector<int> senddata,lengthr,recvdata;
				if(BoundaryRemove[static_cast<size_t>(index)].empty())
					MPI_Send(&trash,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
				{
					vector<int> lengths(BoundaryRemove[static_cast<size_t>(index)].size());
					for(size_t j=0;j<BoundaryNeigh[static_cast<size_t>(index)].size();++j)
					{
					  lengths[j]=static_cast<int>(BoundaryNeigh[static_cast<size_t>(index)][j].size());
						for(int k=0;k<lengths[j];++k)
							senddata.push_back(BoundaryNeigh[static_cast<size_t>(index)][j][k]);
					}
					MPI_Send(&lengths[0],static_cast<int>(lengths.size()),MPI_INT,procorder[i],0,
						MPI_COMM_WORLD);
				}
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&trashr,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					lengthr.resize(static_cast<size_t>(count));
					MPI_Recv(&lengthr[0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if(!BoundaryRemove[static_cast<size_t>(index)].empty())
					MPI_Send(&senddata[0],static_cast<int>(senddata.size()),MPI_INT,procorder[i],0,MPI_COMM_WORLD);
				if(!lengthr.empty())
				{
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					recvdata.resize(static_cast<size_t>(count));
					MPI_Recv(&recvdata[0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				int loc=0;
				// Reorganize the data
				for(size_t j=0;j<lengthr.size();++j)
				{
					vector<int> itemp(recvdata.begin()+loc,recvdata.begin()+loc+
						lengthr[j]);
					loc+=lengthr[j];
					RecvNeigh[static_cast<size_t>(index)].push_back(itemp);
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
					RecvPoints[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
					MPI_Recv(&RecvPoints[static_cast<size_t>(index)][0],count,MPI_INT,procorder[i],0,
						MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if(BoundaryRemove[static_cast<size_t>(index)].empty())
					MPI_Send(&trash,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
				  MPI_Send(&BoundaryRemove[static_cast<size_t>(index)][0],static_cast<int>(BoundaryRemove[static_cast<size_t>(index)].size()),
					MPI_INT,procorder[i],0,MPI_COMM_WORLD);
				MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
				if(stat.MPI_TAG==1)
					MPI_Recv(&trashr,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				else
				{
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					lengthr.resize(static_cast<size_t>(count));
					MPI_Recv(&lengthr[0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if(BoundaryRemove[static_cast<size_t>(index)].empty())
					MPI_Send(&trash,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
				else
				{
					vector<int> lengths(BoundaryRemove[static_cast<size_t>(index)].size());
					for(size_t j=0;j<BoundaryNeigh[static_cast<size_t>(index)].size();++j)
					{
					  lengths[j]=static_cast<int>(BoundaryNeigh[static_cast<size_t>(index)][j].size());
						for(int k=0;k<lengths[j];++k)
							senddata.push_back(BoundaryNeigh[static_cast<size_t>(index)][j][k]);
					}
					MPI_Send(&lengths[0],static_cast<int>(lengths.size()),MPI_INT,procorder[i],0,
						MPI_COMM_WORLD);
				}
				if(!lengthr.empty())
				{
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
					int count;
					MPI_Get_count(&stat,MPI_INT,&count);
					recvdata.resize(static_cast<size_t>(count));
					MPI_Recv(&recvdata[0],count,MPI_INT,procorder[i],0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				if(!BoundaryRemove[static_cast<size_t>(index)].empty())
					MPI_Send(&senddata[0],static_cast<int>(senddata.size()),MPI_INT,procorder[i],0,MPI_COMM_WORLD);
				int loc=0;
				// Reorganize the data
				for(size_t j=0;j<lengthr.size();++j)
				{
					vector<int> itemp(recvdata.begin()+loc,recvdata.begin()+loc+
						lengthr[j]);
					loc+=lengthr[j];
					RecvNeigh[static_cast<size_t>(index)].push_back(itemp);
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
	int nremove=static_cast<int>(toremove.size());
	vector<int> proclist=tess.GetDuplicatedProcs();
	const int rank = get_mpi_rank();
	const int worldsize = get_mpi_size();
	const vector<int> procorder=GetProcOrder(rank,worldsize);
	int nlist=static_cast<int>(proclist.size());
	vector<vector<int> > bremove(static_cast<size_t>(nlist));
	vector<vector<double> > bmerit(static_cast<size_t>(nlist));
	vector<vector<int> > sentpoints=tess.GetDuplicatedPoints();
	vector<vector<int> > sort_indeces(static_cast<size_t>(nlist));
	for(size_t i=0;i<static_cast<size_t>(nlist);++i)
		if(!sentpoints[i].empty())
		{
			sort_index(sentpoints[i],sort_indeces[i]);
			sort(sentpoints[i].begin(),sentpoints[i].end());
		}
	for(size_t i=0;i<static_cast<size_t>(nremove);++i)
		{
		  for(size_t j=0;j<static_cast<size_t>(nlist);++j)
			{
				if(!binary_search(sentpoints[j].begin(),sentpoints[j].end(),
					toremove[i]))
					continue;
				int index=static_cast<int>(lower_bound(sentpoints[j].begin(),sentpoints[j].end(),
								       toremove[i])-sentpoints[j].begin());
				if(index<static_cast<int>(sentpoints[j].size()))
				{
					bremove[j].push_back(sort_indeces[j][static_cast<size_t>(index)]);
					bmerit[j].push_back(merit[i]);
				}
			}
		}

		MPI_Status status;
		// Send/Recv the data
		vector<vector<int> > recvindex(nlist); // the index in the Nghost vector
		vector<vector<double> > recvmerit(static_cast<size_t>(nlist));
		int temp;
		for(size_t i=0;i<procorder.size();++i)
		{
		  int index=static_cast<int>(Find(proclist.begin(),proclist.end(),procorder[i])
					     -proclist.begin());
			if(index<nlist)
			{
				if(rank<procorder[i])
				{
					if(bremove[static_cast<size_t>(index)].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
					  MPI_Send(&bremove[static_cast<size_t>(index)][0],static_cast<int>(bremove[static_cast<size_t>(index)].size()),
						MPI_INT,procorder[i],0,MPI_COMM_WORLD);
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					if(status.MPI_TAG==1)
						MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					else
					{
						int count;
						MPI_Get_count(&status,MPI_INT,&count);
						recvindex[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
						MPI_Recv(&recvindex[static_cast<size_t>(index)][0],count,MPI_INT,procorder[i],0,
							MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
					if(bmerit[static_cast<size_t>(index)].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
					  MPI_Send(&bmerit[static_cast<size_t>(index)][0],static_cast<int>(bmerit[static_cast<size_t>(index)].size()),
						MPI_DOUBLE,procorder[i],0,MPI_COMM_WORLD);
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					if(status.MPI_TAG==1)
						MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					else
					{
						int count;
						MPI_Get_count(&status,MPI_DOUBLE,&count);
						recvmerit[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
						MPI_Recv(&recvmerit[static_cast<size_t>(index)][0],count,MPI_DOUBLE,procorder[i],0,
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
						recvindex[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
						MPI_Recv(&recvindex[static_cast<size_t>(index)][0],count,MPI_INT,procorder[i],0,
							MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
					if(bremove[static_cast<size_t>(index)].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
					  MPI_Send(&bremove[static_cast<size_t>(index)][0],static_cast<int>(bremove[static_cast<size_t>(index)].size()),
						MPI_INT,procorder[i],0,MPI_COMM_WORLD);
					MPI_Probe(procorder[i],MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					if(status.MPI_TAG==1)
						MPI_Recv(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					else
					{
						int count;
						MPI_Get_count(&status,MPI_DOUBLE,&count);
						recvmerit[static_cast<size_t>(index)].resize(static_cast<size_t>(count));
						MPI_Recv(&recvmerit[static_cast<size_t>(index)][0],count,MPI_DOUBLE,procorder[i],0,
							MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					}
					if(bmerit[static_cast<size_t>(index)].empty())
						MPI_Send(&temp,1,MPI_INT,procorder[i],1,MPI_COMM_WORLD);
					else
					  MPI_Send(&bmerit[static_cast<size_t>(index)][0],static_cast<int>(bmerit[static_cast<size_t>(index)].size()),
						MPI_DOUBLE,procorder[i],0,MPI_COMM_WORLD);
				}
			}
		}

		vector<int> res,bad;
		vector<vector<int> > const& Nghost=tess.GetGhostIndeces();
		vector<vector<int> > const& DupPoints=tess.GetDuplicatedPoints();
		for(size_t i=0;i<static_cast<size_t>(nlist);++i)
		{
			if(recvindex[i].empty()||sentpoints[i].empty())
				continue;
			vector<int> Nghostindex(recvindex[i].size());
			for(size_t j=0;j<Nghostindex.size();++j)
				Nghostindex[j]=Nghost[i][recvindex[i][j]];
			vector<int> indeces;
			sort_index(Nghostindex,indeces);
			sort(Nghostindex.begin(),Nghostindex.end());
			for(size_t j=0;j<bremove[i].size();++j)
			{
				vector<int> neigh=tess.GetNeighbors(DupPoints[i][bremove[i][j]]);
				for(int k=0;k<static_cast<int>(neigh.size());++k)
				{
					if(!binary_search(Nghostindex.begin(),Nghostindex.end(),neigh[k]))
						continue;
					int index=static_cast<int>(lower_bound(Nghostindex.begin(),Nghostindex.end(),neigh[k])-
								   Nghostindex.begin());
					if(index<static_cast<int>(Nghostindex.size()))
						if(recvmerit[i][indeces[static_cast<size_t>(index)]]>=bmerit[i][j])
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

namespace {
  size_t sum_sizes(const vector<vector<int> >& vvi)
  {
    size_t sum = 0;
    for(size_t i=0;i<vvi.size();++i)
      sum += vvi[i].size();
    return sum;
  }
}

void PeriodicUpdateCells(vector<Primitive> &cells,vector<vector<double> > &tracers,
	vector<size_t> &customevolutions,vector<vector<int> > const& sentcells,
	int npoints)
{
	if(sentcells.empty())
		return;
	const bool traceractive = !tracers.empty() && !tracers.at(0).empty();
	const int totalsent = static_cast<int>(sum_sizes(sentcells));
	cells.resize(static_cast<size_t>(npoints-totalsent));
	customevolutions.resize(static_cast<size_t>(npoints-totalsent));
	if(traceractive)
	  tracers.resize(static_cast<size_t>(npoints-totalsent));
	for(size_t i=0;i<sentcells.size();++i)
	{
		if(sentcells[i].empty())
			continue;
		insert_all_to_back(cells,
				   VectorValues(cells,sentcells[i]));
		insert_all_to_back(customevolutions,
				   VectorValues(customevolutions,sentcells[i]));
		if(traceractive)
		  insert_all_to_back(tracers,
				     VectorValues(tracers,sentcells[i]));
	}
}

void PeriodicVelocityExchange(vector<Vector2D> &vel,
	vector<vector<int> > const& sentcells,int npoints)
{
  const int totalsent = static_cast<int>(sum_sizes(sentcells));
	vel.resize(static_cast<size_t>(npoints-totalsent));

	for(size_t i=0;i<sentcells.size();++i)
	{
	  if(!sentcells[i].empty())
	    insert_all_to_back(vel,
			       VectorValues(vel,sentcells[i]));
	}
}

void PeriodicGradExchange(vector<ReducedPrimitiveGradient2D> &grad,
	vector<vector<int> > const& sentcells,int npoints)
{
  const int totalsent = static_cast<int>(sum_sizes(sentcells));
  grad.resize(static_cast<size_t>(npoints-totalsent));

	for(size_t i=0;i<sentcells.size();++i)
	{
	  if(!sentcells[i].empty())
	    insert_all_to_back(grad,VectorValues(grad,sentcells[i]));
	}
}
